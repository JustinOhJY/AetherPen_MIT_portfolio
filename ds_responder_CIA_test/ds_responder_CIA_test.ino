#include "dw3000.h"
 
#define APP_NAME "DS TWR RESP v1.0"
 
// connection pins
// connection pins
const uint8_t PIN_RST = 15; // reset pin
const uint8_t PIN_IRQ = 4; // irq pin
const uint8_t PIN_SS = 5; // spi select pin
 
/* Default communication configuration. We use default non-STS DW mode. */
static dwt_config_t config = {
  5,               /* Channel number. */
  DWT_PLEN_128,    /* Preamble length. Used in TX only. */
  DWT_PAC8,        /* Preamble acquisition chunk size. Used in RX only. */
  9,               /* TX preamble code. Used in TX only. */
  9,               /* RX preamble code. Used in RX only. */
  1,               /* 0 to use standard 8 symbol SFD, 1 to use non-standard 8 symbol, 2 for non-standard 16 symbol SFD and 3 for 4z 8 symbol SDF type */
  DWT_BR_6M8,      /* Data rate. */
  DWT_PHRMODE_STD, /* PHY header mode. */
  DWT_PHRRATE_STD, /* PHY header rate. */
  (129 + 8 - 8),   /* SFD timeout (preamble length + 1 + SFD length - PAC size). Used in RX only. */
  DWT_STS_MODE_OFF, /* STS disabled */
  DWT_STS_LEN_64,/* STS length see allowed values in Enum dwt_sts_lengths_e */
  DWT_PDOA_M0      /* PDOA mode off */
};

#define IP_CONF_REG  0x0E   // Register file for CIA registers
#define IP_CONF_SUB  0x0C   // Sub-register offset for IP_CONF

// Field definitions:
// IP_NTM occupies bits 0-4 (5 bits)
#define IP_NTM_MASK  (0x1F << 0)   // 0x1F = 0b11111
#define IP_NTM_SHIFT 0

// IP_PMULT occupies bits 5-6 (2 bits)
#define IP_PMULT_MASK  (0x3 << 5)   // 0x3 = 0b11; shifted left 5 gives 0x60 (bits 5-6)
#define IP_PMULT_SHIFT 5

// IP_RTM occupies bits 16-20 (5 bits)
#define IP_RTM_MASK  (0x1F << 16)   // 0x1F = 0b11111; shifted left 16 places
#define IP_RTM_SHIFT 16
/* Inter-ranging delay period, in milliseconds. */
#define RNG_DELAY_MS 500

/* Default antenna delay values for 64 MHz PRF. See NOTE 1 below. */
#define TX_ANT_DLY 16385
#define RX_ANT_DLY 16385

/* Frames used in the ranging process. See NOTE 2 below. */
static uint8_t rx_poll_msg[] = {0x41, 0x88, 0, 0xCA, 0xDE, 'W', 'A', 'V', 'E', 0x21, 0, 0};
static uint8_t tx_resp_msg[] = {0x41, 0x88, 0, 0xCA, 0xDE, 'V', 'E', 'W', 'A', 0x10, 0x02, 0, 0, 0, 0};
static uint8_t rx_final_msg[] = {0x41, 0x88, 0, 0xCA, 0xDE, 'W', 'A', 'V', 'E', 0x23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/* Length of the common part of the message (up to and including the
 * function code, see NOTE 2 below). */
#define ALL_MSG_COMMON_LEN 10

/* Index to access some of the fields in the frames involved in the process. */
#define ALL_MSG_SN_IDX 2
#define FINAL_MSG_POLL_TX_TS_IDX 10
#define FINAL_MSG_RESP_RX_TS_IDX 14
#define FINAL_MSG_FINAL_TX_TS_IDX 18

/* Frame sequence number, incremented after each transmission. */
static uint8_t frame_seq_nb = 0;

/* Buffer to store received messages.
 * Its size is adjusted to longest frame that this example code is
 * supposed to handle. */
#define RX_BUF_LEN 24
static uint8_t rx_buffer[RX_BUF_LEN];

/* Hold copy of status register state here for reference so that it can be
 * examined at a debug breakpoint. */
static uint32_t status_reg = 0;

/* Delay between frames, in UWB microseconds. See NOTE 4 below. */
/* This is the delay from Frame RX timestamp to TX reply timestamp used for
 * calculating/setting the DW IC's delayed TX function. This includes the
 * frame length of approximately 190 us with above configuration. */
#define POLL_RX_TO_RESP_TX_DLY_UUS 900

/* This is the delay from the end of the frame transmission to the enable of
 * the receiver, as programmed for the DW IC's wait for response feature. */
#define RESP_TX_TO_FINAL_RX_DLY_UUS 500

/* Receive final timeout. See NOTE 5 below. */
#define FINAL_RX_TIMEOUT_UUS 220

/* Preamble timeout, in multiple of PAC size. See NOTE 6 below. */
#define PRE_TIMEOUT 5

/* Timestamps of frames transmission/reception. */
static uint64_t poll_rx_ts;
static uint64_t resp_tx_ts;
static uint64_t final_rx_ts;

/* Hold copies of computed time of flight and distance here for reference
 * so that it can be examined at a debug breakpoint. */
static double tof;
static double distance;

/* Values for the PG_DELAY and TX_POWER registers reflect the bandwidth and
 * power of the spectrum at the current temperature. These values can be
 * calibrated prior to taking reference measurements.
 * See NOTE 2 below. */
extern dwt_txconfig_t txconfig_options;

/*! ---------------------------------------------------------------------------
 * @fn ds_twr_responder()
 *
 * @brief Application entry point.
 *
 * @param  none
 *
 * @return none
 */
void setup() 
{
  UART_init();
  test_run_info((unsigned char *)APP_NAME);
  spiBegin(PIN_IRQ, PIN_RST);
  spiSelect(PIN_SS);
    /* Display application name on LCD. */
   // LOG_INF(APP_NAME);

    /* Configure SPI rate, DW3000 supports up to 38 MHz */
   // port_set_dw_ic_spi_fastrate();

    /* Reset DW IC */
   // reset_DWIC();

    Sleep(2);

    /* Need to make sure DW IC is in IDLE_RC before proceeding */
    while (!dwt_checkidlerc()) {
          UART_puts("IDLE FAILED\r\n");
    while (1) ;

     };

    if (dwt_initialise(DWT_DW_INIT) == DWT_ERROR) {
       // LOG_ERR("INIT FAILED");
        while (1) { /* spin */ };
    }

    /* Configure DW IC. See NOTE 15 below. */
    if (dwt_configure(&config)) {
       // LOG_ERR("CONFIG FAILED");
        while (1) { /* spin */ };
    }

// *** ADDED: Definitions for IP_CONF, IP_NTM, IP_PMULT, and IP_RTM ***
#define IP_CONF_REG  0x0E   // Register file for CIA parameters
#define IP_CONF_SUB  0x0C   // Sub-register offset for IP_CONF

// According to your updated info:
// IP_NTM occupies bits 0-4 (5 bits): mask = 0x1F (0b11111)
// IP_PMULT occupies bits 5-6 (2 bits): mask = 0x3 << 5 = 0x60
// IP_RTM occupies bits 16-20 (5 bits): mask = 0x1F << 16
#define IP_NTM_MASK  (0x1F << 0)
#define IP_NTM_SHIFT 0

#define IP_PMULT_MASK  (0x3 << 5)
#define IP_PMULT_SHIFT 5

#define IP_RTM_MASK  (0x1F << 16)
#define IP_RTM_SHIFT 16

// *** ADDED: In setup(), after dwt_configure(&config) succeeds, update thresholds:
{
    // Read the current 32-bit IP_CONF value from file 0x0E, sub-register 0x0C.
    uint32_t ip_conf_val = dwt_read32bitoffsetreg(IP_CONF_REG, IP_CONF_SUB);

    // Clear existing bits for IP_NTM, IP_PMULT, and IP_RTM.
    ip_conf_val &= ~(IP_NTM_MASK | IP_PMULT_MASK | IP_RTM_MASK);

    // Set new values:
    // Set IP_NTM to 0x1C (28 decimal) to raise the noise floor multiplier.
    ip_conf_val |= ((0x10 << IP_NTM_SHIFT) & IP_NTM_MASK);

    // Set IP_PMULT to 0x3 (max for 2 bits).
    ip_conf_val |= ((0x2 << IP_PMULT_SHIFT) & IP_PMULT_MASK);

    // Set IP_RTM to 0x18 (24 decimal) to adjust the search window.
    ip_conf_val |= ((0x12 << IP_RTM_SHIFT) & IP_RTM_MASK);

    // Write the updated value back to the IP_CONF register.
    dwt_write32bitoffsetreg(IP_CONF_REG, IP_CONF_SUB, ip_conf_val);

    UART_puts("CIA thresholds updated: IP_NTM=0x1C, IP_PMULT=0x3, IP_RTM=0x18\r\n");
}

    /* Configure the TX spectrum parameters (power, PG delay and PG count) */
    dwt_configuretxrf(&txconfig_options);

    /* Apply default antenna delay value. See NOTE 1 below. */
    dwt_setrxantennadelay(RX_ANT_DLY);
    dwt_settxantennadelay(TX_ANT_DLY);

    /* Next can enable TX/RX states output on GPIOs 5 and 6 to help debug,
     * and also TX/RX LEDs.
     * Note, in real low power applications the LEDs should not be used.
     */
    dwt_setlnapamode(DWT_LNA_ENABLE | DWT_PA_ENABLE);

    dwt_setleds(DWT_LEDS_ENABLE | DWT_LEDS_INIT_BLINK);
}    

    /* Loop forever responding to ranging requests. */
    void loop () 
    {while (1){

        dwt_setpreambledetecttimeout(0);
        /* Clear reception timeout to start next ranging process. */
        dwt_setrxtimeout(0);

        /* Activate reception immediately. */
        dwt_rxenable(DWT_START_RX_IMMEDIATE);

        /* Poll for reception of a frame or error/timeout. See NOTE 8 below. */
        while (!((status_reg = dwt_read32bitreg(SYS_STATUS_ID)) & (SYS_STATUS_RXFCG_BIT_MASK |SYS_STATUS_ALL_RX_TO |SYS_STATUS_ALL_RX_ERR)))
        { /* spin */ };

        if (status_reg & SYS_STATUS_RXFCG_BIT_MASK) {

            /* Clear good RX frame event in the DW IC status register. */
            dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_RXFCG_BIT_MASK);

            /* A frame has been received, read it into the local buffer. */
            uint32_t frame_len = dwt_read32bitreg(RX_FINFO_ID) & FRAME_LEN_MAX_EX;
            if (frame_len <= RX_BUF_LEN) {
                dwt_readrxdata(rx_buffer, frame_len, 0);
            }

            /* Check that the frame is a poll sent by "DS TWR initiator" example.
             * As the sequence number field of the frame is not relevant, it
             * is cleared to simplify the validation of the frame.
             */
            rx_buffer[ALL_MSG_SN_IDX] = 0;
            if (memcmp(rx_buffer, rx_poll_msg, ALL_MSG_COMMON_LEN) == 0) {

                uint32_t resp_tx_time;

                /* Retrieve poll reception timestamp. */
                poll_rx_ts = get_rx_timestamp_u64();

                /* Set send time for response. See NOTE 9 below. */
                resp_tx_time = (poll_rx_ts + (POLL_RX_TO_RESP_TX_DLY_UUS * UUS_TO_DWT_TIME)) >> 8;
                dwt_setdelayedtrxtime(resp_tx_time);

                /* Set expected delay and timeout for final message reception. See NOTE 4 and 5 below. */
                dwt_setrxaftertxdelay(RESP_TX_TO_FINAL_RX_DLY_UUS);
                dwt_setrxtimeout(FINAL_RX_TIMEOUT_UUS);

                /* Set preamble timeout for expected frames. See NOTE 6 below. */
                dwt_setpreambledetecttimeout(PRE_TIMEOUT);

                /* Write and send the response message. See NOTE 10 below.*/
                tx_resp_msg[ALL_MSG_SN_IDX] = frame_seq_nb;
                dwt_writetxdata(sizeof(tx_resp_msg), tx_resp_msg, 0); /* Zero offset in TX buffer. */
                dwt_writetxfctrl(sizeof(tx_resp_msg), 0, 1); /* Zero offset in TX buffer, ranging. */
                int ret = dwt_starttx(DWT_START_TX_DELAYED | DWT_RESPONSE_EXPECTED);

                /* If dwt_starttx() returns an error, abandon this ranging
                 * exchange and proceed to the next one. See NOTE 11 below. */
                if (ret == DWT_ERROR) {
                    continue;
                }

                /* Poll for reception of expected "final" frame or error/timeout.
                 * See NOTE 8 below.
                 */
                while (!((status_reg = dwt_read32bitreg(SYS_STATUS_ID)) &
                                                       (SYS_STATUS_RXFCG_BIT_MASK |
                                                        SYS_STATUS_ALL_RX_TO |
                                                        SYS_STATUS_ALL_RX_ERR)))
                { /* spin */ };

                /* Increment frame sequence number after transmission of the
                 * response message (modulo 256).
                 */
                frame_seq_nb++;

                if (status_reg & SYS_STATUS_RXFCG_BIT_MASK) {

                    /* Clear good RX frame event and TX frame sent in the DW IC status register. */
                    dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_RXFCG_BIT_MASK | SYS_STATUS_TXFRS_BIT_MASK);

                    /* A frame has been received, read it into the local buffer. */
                    frame_len = dwt_read32bitreg(RX_FINFO_ID) & FRAME_LEN_MAX_EX;
                    if (frame_len <= RX_BUF_LEN) {
                        dwt_readrxdata(rx_buffer, frame_len, 0);
                    }

                    /* Check that the frame is a final message sent by
                     * "DS TWR initiator" example.
                     * As the sequence number field of the frame is not used in
                     * this example, it can be zeroed to ease the validation of
                     * the frame.
                     */
                    rx_buffer[ALL_MSG_SN_IDX] = 0;
                    if (memcmp(rx_buffer, rx_final_msg, ALL_MSG_COMMON_LEN) == 0) {

                        uint32_t poll_tx_ts, resp_rx_ts, final_tx_ts;
                        uint32_t poll_rx_ts_32, resp_tx_ts_32, final_rx_ts_32;
                        double Ra, Rb, Da, Db;
                        int64_t tof_dtu;

                        /* Retrieve response transmission and final
                         * reception timestamps. */
                        resp_tx_ts = get_tx_timestamp_u64();
                        final_rx_ts = get_rx_timestamp_u64();

                        /* Get timestamps embedded in the final message. */
                        final_msg_get_ts(&rx_buffer[FINAL_MSG_POLL_TX_TS_IDX], &poll_tx_ts);
                        final_msg_get_ts(&rx_buffer[FINAL_MSG_RESP_RX_TS_IDX], &resp_rx_ts);
                        final_msg_get_ts(&rx_buffer[FINAL_MSG_FINAL_TX_TS_IDX], &final_tx_ts);

                        /* Compute time of flight. 32-bit subtractions give
                         * correct answers even if clock has wrapped.
                         * See NOTE 12 below. */
                        poll_rx_ts_32 = (uint32_t)poll_rx_ts;
                        resp_tx_ts_32 = (uint32_t)resp_tx_ts;
                        final_rx_ts_32 = (uint32_t)final_rx_ts;
                        Ra = (double)(resp_rx_ts - poll_tx_ts);
                        Rb = (double)(final_rx_ts_32 - resp_tx_ts_32);
                        Da = (double)(final_tx_ts - resp_rx_ts);
                        Db = (double)(resp_tx_ts_32 - poll_rx_ts_32);
                        tof_dtu = (int64_t)((Ra * Rb - Da * Db) / (Ra + Rb + Da + Db));

                        tof = tof_dtu * DWT_TIME_UNITS;
                        distance = tof * SPEED_OF_LIGHT;

                        /* Display computed distance. */
                        //char dist[20] = {0};
                    // snprintf(dist_str, sizeof(dist_str), " %3.1f ", distance);
                    // test_run_info((unsigned char *)dist_str);
                    Serial.println("0 "+String(distance*100)+" 200"); // LowerLimit Reading UpperLimit
                        //LOG_INF("%s", log_strdup(dist));

                        /* As DS-TWR initiator is waiting for RNG_DELAY_MS
                         * before next poll transmission we can add a delay
                         * here before RX is re-enabled again.
                         */
                        Sleep(RNG_DELAY_MS - 10);  // start couple of ms earlier
                    }
                }
                else {
                    /* Clear RX error/timeout events in the DW IC
                     * status register. */
                    dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR);
                }
                
            }
        }
        else 
        {
            /* Clear RX error/timeout events in the DW IC status register. */
            dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR);
        }
    }
    }