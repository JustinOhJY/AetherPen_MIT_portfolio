#include "dw3000.h"

#define APP_NAME "DS TWR RESP (No Calculation) v1.0"

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
 * The poll message from the tag is 12 bytes (header + 2 extra bytes).
 * The response message is extended to 28 bytes:
 *  - First 10 bytes: common header.
 *  - Next 2 bytes: (if needed)
 *  - Then 8 bytes for T₂ (poll reception timestamp at anchor)
 *  - Then 8 bytes for T₃ (response transmission timestamp at anchor)
 */
static uint8_t rx_poll_msg[] = {0x41, 0x88, 0, 0xCA, 0xDE, 'W','A','V','E', 0x21, 0, 0};

#define RESP_MSG_LEN 28
static uint8_t tx_resp_msg[RESP_MSG_LEN] = {
  0x41, 0x88, 0, 0xCA, 0xDE, 'V','E','W','A', 0x10, 0x02, 0,
  // Placeholders for T₂ (8 bytes)
  0,0,0,0,0,0,0,0,
  // Placeholders for T₃ (8 bytes)
  0,0,0,0,0,0,0,0
};

/* Length of the common part of the message. */
#define ALL_MSG_COMMON_LEN 10

/* Index definitions */
#define ALL_MSG_SN_IDX 2
// Starting indexes for the timestamps in the response message:
#define RESP_T2_IDX 12
#define RESP_T3_IDX 20

/* Buffer to store received messages. */
#define RX_BUF_LEN  32
static uint8_t rx_buffer[RX_BUF_LEN];

/* Hold copy of status register state. */
static uint32_t status_reg = 0;

/* Timing delay for response transmission in UWB microseconds. */
#define POLL_RX_TO_RESP_TX_DLY_UUS 900
#define PRE_TIMEOUT                5

/* Timestamp for poll reception at anchor (T₂). */
static uint64_t poll_rx_ts;

/* External TX configuration parameters. */
extern dwt_txconfig_t txconfig_options;

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
    dwt_setpreambledetecttimeout(PRE_TIMEOUT);
    dwt_setlnapamode(DWT_LNA_ENABLE | DWT_PA_ENABLE);
    dwt_setleds(DWT_LEDS_ENABLE | DWT_LEDS_INIT_BLINK);
}

void loop() {
    while (1) {
        // Prepare the receiver to listen for a poll message.
        dwt_setpreambledetecttimeout(0);
        dwt_setrxtimeout(0);
        dwt_rxenable(DWT_START_RX_IMMEDIATE);
        
        // Wait until a poll is received or a timeout/error occurs.
        while (!( (status_reg = dwt_read32bitreg(SYS_STATUS_ID)) &
                  (SYS_STATUS_RXFCG_BIT_MASK | SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR) ))
        { /* spin */ };

        if (status_reg & SYS_STATUS_RXFCG_BIT_MASK) {
            // Clear the RX good frame event.
            dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_RXFCG_BIT_MASK);

            uint32_t frame_len = dwt_read32bitreg(RX_FINFO_ID) & FRAME_LEN_MAX_EX;
            if (frame_len <= RX_BUF_LEN) {
                dwt_readrxdata(rx_buffer, frame_len, 0);
            }
            // Validate the poll message header (ignore sequence number)
            rx_buffer[ALL_MSG_SN_IDX] = 0;
            if (memcmp(rx_buffer, rx_poll_msg, ALL_MSG_COMMON_LEN) == 0) {
                // Record the poll reception timestamp (T₂).
                poll_rx_ts = get_rx_timestamp_u64();
                
                // Set up delayed transmission for the response.
                uint32_t resp_tx_time = (poll_rx_ts + (POLL_RX_TO_RESP_TX_DLY_UUS * UUS_TO_DWT_TIME)) >> 8;
                dwt_setdelayedtrxtime(resp_tx_time);
                
                // Pre-calculate the anchor's response TX timestamp (T₃).
                uint64_t anchor_resp_tx_ts = (((uint64_t)(resp_tx_time & 0xFFFFFFFEUL)) << 8) + TX_ANT_DLY;
                
                // Update sequence number if desired.
                tx_resp_msg[ALL_MSG_SN_IDX] = rx_poll_msg[ALL_MSG_SN_IDX];
                
                // Write T₂ and T₃ into the response message.
                final_msg_set_ts(&tx_resp_msg[RESP_T2_IDX], poll_rx_ts);
                final_msg_set_ts(&tx_resp_msg[RESP_T3_IDX], anchor_resp_tx_ts);
                
                // Write and send the response message.
                dwt_writetxdata(sizeof(tx_resp_msg), tx_resp_msg, 0);
                dwt_writetxfctrl(sizeof(tx_resp_msg) + FCS_LEN, 0, 1);
                
                int ret = dwt_starttx(DWT_START_TX_DELAYED);
                if (ret == DWT_ERROR) {
                    continue; // Abandon this exchange if TX fails.
                }
                // Wait for the TX frame sent event.
                while (!(dwt_read32bitreg(SYS_STATUS_ID) & SYS_STATUS_TXFRS_BIT_MASK))
                { /* spin */ };
                dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_TXFRS_BIT_MASK);
            }
        } else {
            // Clear any RX error/timeout events.
            dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR);
        }
    }
}
