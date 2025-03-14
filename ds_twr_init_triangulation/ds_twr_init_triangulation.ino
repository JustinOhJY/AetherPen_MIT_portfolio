#include "dw3000.h"

#define APP_NAME "DS TWR INIT (Distance on Tag) v1.0"

// Connection pins
const uint8_t PIN_RST = 15; // reset pin
const uint8_t PIN_IRQ = 4; // irq pin
const uint8_t PIN_SS  = 5;  // spi select pin

/* Default communication configuration. We use default non-STS DW mode. */
static dwt_config_t config = {
    5,               /* Channel number. */
    DWT_PLEN_128,    /* Preamble length. */
    DWT_PAC8,        /* Preamble acquisition chunk size. */
    9,               /* TX preamble code. */
    9,               /* RX preamble code. */
    1,               /* SFD mode */
    DWT_BR_6M8,      /* Data rate. */
    DWT_PHRMODE_STD, /* PHY header mode. */
    DWT_PHRRATE_STD, /* PHY header rate. */
    (129 + 8 - 8),   /* SFD timeout */
    DWT_STS_MODE_OFF,/* STS disabled */
    DWT_STS_LEN_64,  /* STS length */
    DWT_PDOA_M0      /* PDOA mode off */
};

/* Inter-ranging delay period, in milliseconds. */
#define RNG_DELAY_MS 500

/* Default antenna delay values for 64 MHz PRF. */
#define TX_ANT_DLY 16385
#define RX_ANT_DLY 16385

/* Frames used in the ranging process.  
 * The poll message is 10 bytes. The response message is expected to be 28 bytes:
 *  - First 10 bytes: header (common part)
 *  - Next 2 bytes: extra header fields (if any)
 *  - Then 8 bytes for T₂ (anchor poll-reception timestamp)
 *  - Then 8 bytes for T₃ (anchor response-transmission timestamp)
 */
static uint8_t tx_poll_msg[] = {0x41, 0x88, 0, 0xCA, 0xDE, 'W','A','V','E', 0x21};

#define RESP_MSG_LEN 28
static uint8_t rx_resp_msg[RESP_MSG_LEN] = {
  0x41, 0x88, 0, 0xCA, 0xDE, 'V','E','W','A', 0x10, 0x02, 0,
  // 16 bytes for anchor's T2 and T3 (initialized to zeros)
  0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0
};

static uint8_t frame_seq_nb = 0;


/* Length of the common part of the message. */
#define ALL_MSG_COMMON_LEN 10

/* Index definitions */
#define ALL_MSG_SN_IDX 2
// Starting indexes of the timestamps in the response message:
#define RESP_T2_IDX 12
#define RESP_T3_IDX 20

/* Buffer to store the received response message. */
#define RX_BUF_LEN  RESP_MSG_LEN
static uint8_t rx_buffer[RX_BUF_LEN];

/* Hold copy of status register state for reference. */
static uint32_t status_reg = 0;

/* Timing delays in UWB microseconds. */
#define POLL_TX_TO_RESP_RX_DLY_UUS 700
#define RESP_RX_TIMEOUT_UUS        300
#define PRE_TIMEOUT                5

/* Timestamps (in device time units). */
static uint64_t poll_tx_ts;
static uint64_t resp_rx_ts;

// Speed of light (m/s) for converting time-of-flight to distance.
#define SPEED_OF_LIGHT 299702547.0

/* External TX configuration parameters. */
extern dwt_txconfig_t txconfig_options;

/* Helper function to extract a 64-bit timestamp from the given buffer.
 * (Assumes little-endian format.)
 */
static uint64_t get_timestamp_from_buffer(uint8_t *buffer) {
    uint64_t ts = 0;
    for (int i = 7; i >= 0; i--) {
        ts = (ts << 8) | buffer[i];
    }
    return ts;
}

/*! ---------------------------------------------------------------------------
 * @fn ds_twr_initiator()
 *
 * @brief Modified DS-TWR initiator with distance calculation on tag.
 *
 * @param  none
 *
 * @return none
 */
void setup() {
    UART_init();
    test_run_info((unsigned char *)APP_NAME);
    spiBegin(PIN_IRQ, PIN_RST);
    spiSelect(PIN_SS);

    Sleep(2);
    while (!dwt_checkidlerc()) {
        UART_puts("IDLE FAILED\r\n");
        while (1) ;
    }
    if (dwt_initialise(DWT_DW_INIT) == DWT_ERROR) {
        UART_puts("INIT FAILED\r\n");
        while (1) { }
    }
    if (dwt_configure(&config)) {
        UART_puts("CONFIG FAILED\r\n");
        while (1) { }
    }
    dwt_configuretxrf(&txconfig_options);
    dwt_setrxantennadelay(RX_ANT_DLY);
    dwt_settxantennadelay(TX_ANT_DLY);
    dwt_setrxaftertxdelay(POLL_TX_TO_RESP_RX_DLY_UUS);
    dwt_setrxtimeout(RESP_RX_TIMEOUT_UUS);
    dwt_setpreambledetecttimeout(PRE_TIMEOUT);
    dwt_setlnapamode(DWT_LNA_ENABLE | DWT_PA_ENABLE);
    dwt_setleds(DWT_LEDS_ENABLE | DWT_LEDS_INIT_BLINK);
}

void loop() {
    while (1) {
        // Prepare poll message with updated sequence number
        tx_poll_msg[ALL_MSG_SN_IDX] = frame_seq_nb;
        dwt_writetxdata(sizeof(tx_poll_msg), tx_poll_msg, 0);
        dwt_writetxfctrl(sizeof(tx_poll_msg) + FCS_LEN, 0, 1);

        // Start transmission (T1) and expect a response.
        dwt_starttx(DWT_START_TX_IMMEDIATE | DWT_RESPONSE_EXPECTED);
        poll_tx_ts = get_tx_timestamp_u64();

        // Wait for the response (which now includes T2 and T3).
        while (!( (status_reg = dwt_read32bitreg(SYS_STATUS_ID)) &
                  (SYS_STATUS_RXFCG_BIT_MASK | SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR) ))
        { /* spin */ };

        frame_seq_nb++;  // Increment sequence number

        if (status_reg & SYS_STATUS_RXFCG_BIT_MASK) {
            dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_RXFCG_BIT_MASK | SYS_STATUS_TXFRS_BIT_MASK);
            uint32_t frame_len = dwt_read32bitreg(RX_FINFO_ID) & FRAME_LEN_MAX_EX;
            if (frame_len <= RX_BUF_LEN) {
                dwt_readrxdata(rx_buffer, frame_len, 0);
            }
            // Validate response header (ignore sequence number)
            rx_buffer[ALL_MSG_SN_IDX] = 0;
            if (memcmp(rx_buffer, rx_resp_msg, ALL_MSG_COMMON_LEN) == 0) {
                // Record response reception timestamp (T4)
                resp_rx_ts = get_rx_timestamp_u64();

                // Extract anchor's timestamps T2 and T3 from response message.
                uint64_t anchor_poll_rx_ts = get_timestamp_from_buffer(&rx_buffer[RESP_T2_IDX]);
                uint64_t anchor_resp_tx_ts = get_timestamp_from_buffer(&rx_buffer[RESP_T3_IDX]);

                // Compute time-of-flight:
                // ToF = ((T4 - T1) - (T3 - T2)) / 2
                int64_t round_trip = (int64_t)(resp_rx_ts - poll_tx_ts);
                int64_t anchor_delay = (int64_t)(anchor_resp_tx_ts - anchor_poll_rx_ts);
                int64_t tof_dtu = (round_trip - anchor_delay) / 2;
                
                // Convert device time units to seconds then compute distance.
                double tof = tof_dtu * DWT_TIME_UNITS;
                double distance = tof * SPEED_OF_LIGHT;
                
                // Output the computed distance (in centimeters).
                Serial.println("Distance (cm): " + String(distance * 100));
            }
        } else {
            // Clear any RX error/timeout events.
            dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR | SYS_STATUS_TXFRS_BIT_MASK);
        }
        
        Sleep(RNG_DELAY_MS);
    }
}
