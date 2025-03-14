#include "dw3000.h"

#define APP_NAME "DS TWR RESP v1.0 - Modified for STS"

// Connection pins
const uint8_t PIN_RST = 15; // Reset pin
const uint8_t PIN_IRQ = 4;  // IRQ pin
const uint8_t PIN_SS  = 5;  // SPI select pin

/* Default communication configuration.
 * Previously, STS was disabled; now we enable it by using DWT_STS_MODE_AUTO.
 */
static dwt_config_t config = {
  5,               /* Channel number. */
  DWT_PLEN_128,    /* Preamble length. Used in TX only. */
  DWT_PAC8,        /* Preamble acquisition chunk size. Used in RX only. */
  9,               /* TX preamble code. */
  9,               /* RX preamble code. */
  1,               /* SFD type (non-standard 8-symbol SFD). */
  DWT_BR_6M8,      /* Data rate. */
  DWT_PHRMODE_STD, /* PHY header mode. */
  DWT_PHRRATE_STD, /* PHY header rate. */
  (129 + 8 - 8),   /* SFD timeout. */
  DWT_STS_MODE_1, /* Enable STS instead of off */
  DWT_STS_LEN_64,  /* STS length. */
  DWT_PDOA_M0      /* PDOA mode off. */
};

/* Inter-ranging delay period, in milliseconds. */
#define RNG_DELAY_MS 500

/* Default antenna delay values for 64 MHz PRF. */
#define TX_ANT_DLY 16385
#define RX_ANT_DLY 16385

/* Frames used in the ranging process. */
static uint8_t rx_poll_msg[] = { 0x41, 0x88, 0, 0xCA, 0xDE, 'W', 'A', 'V', 'E', 0x21, 0, 0 };
static uint8_t tx_resp_msg[] = { 0x41, 0x88, 0, 0xCA, 0xDE, 'V', 'E', 'W', 'A', 0x10, 0x02, 0, 0, 0, 0 };
static uint8_t rx_final_msg[] = { 0x41, 0x88, 0, 0xCA, 0xDE, 'W', 'A', 'V', 'E', 0x23,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

/* Length of the common part of the message (first 10 bytes). */
#define ALL_MSG_COMMON_LEN 10

/* Index to access fields in the messages. */
#define ALL_MSG_SN_IDX             2
#define FINAL_MSG_POLL_TX_TS_IDX  10
#define FINAL_MSG_RESP_RX_TS_IDX  14
#define FINAL_MSG_FINAL_TX_TS_IDX 18

/* Frame sequence number, incremented after each transmission. */
static uint8_t frame_seq_nb = 0;

/* Buffer to store received messages. */
#define RX_BUF_LEN 24
static uint8_t rx_buffer[RX_BUF_LEN];

/* Copy of status register for debugging. */
static uint32_t status_reg = 0;

/* Timing parameters (in UWB microseconds). */
#define POLL_RX_TO_RESP_TX_DLY_UUS 900
#define RESP_TX_TO_FINAL_RX_DLY_UUS 500
#define FINAL_RX_TIMEOUT_UUS        220
#define PRE_TIMEOUT                 5

/* Timestamps of frame events (device time units). */
static uint64_t poll_rx_ts;
static uint64_t resp_tx_ts;
static uint64_t final_rx_ts;

/* External TX configuration (calibrated per device). */
extern dwt_txconfig_t txconfig_options;

/*
 * --- STS Key and IV Definitions ---
 * These must be identical on both the initiator and responder.
 * The defaults here are as specified in the IEEE 802.15.4z annex.
 */
static dwt_sts_cp_key_t cp_key = {
    0x14EB220F, 0xF86050A8, 0xD1D336AA, 0x14148674
};

static dwt_sts_cp_iv_t cp_iv = {
    0x1F9A3DE4, 0xD37EC3CA, 0xC44FA8FB, 0x362EEB34
};

/* Flag to ensure STS key/IV are loaded only once at startup. */
static uint8_t firstLoopFlag = 0;

/*---------------------------------------------------------------------------
 * DS TWR RESP Setup
 *---------------------------------------------------------------------------
 */
void setup() 
{
  UART_init();
  test_run_info((unsigned char *)APP_NAME);
  spiBegin(PIN_IRQ, PIN_RST);
  spiSelect(PIN_SS);
  
  /* Allow time for DW3000 startup (transition from INIT_RC to IDLE_RC). */
  Sleep(2);
  
  /* Ensure DW3000 is in IDLE state. */
  while (!dwt_checkidlerc())
  {
    UART_puts("IDLE FAILED\r\n");
    while (1) ;
  }
  
  if (dwt_initialise(DWT_DW_INIT) == DWT_ERROR)
  {
    UART_puts("INIT FAILED\r\n");
    while (1) ;
  }
  
  /* Optionally enable LEDs for debug. */
  dwt_setleds(DWT_LEDS_ENABLE | DWT_LEDS_INIT_BLINK);
  
  /* Configure DW3000 with our STS-enabled DS configuration. */
  if (dwt_configure(&config))
  {
    UART_puts("CONFIG FAILED\r\n");
    while (1) ;
  }
  
  /* Configure TX spectrum parameters. */
  dwt_configuretxrf(&txconfig_options);
  
  /* Set default antenna delays. */
  dwt_setrxantennadelay(RX_ANT_DLY);
  dwt_settxantennadelay(TX_ANT_DLY);
  
  /* Enable LNA/PA and optionally debug LEDs. */
  dwt_setlnapamode(DWT_LNA_ENABLE | DWT_PA_ENABLE);
  dwt_setleds(DWT_LEDS_ENABLE | DWT_LEDS_INIT_BLINK);
}

/*---------------------------------------------------------------------------
 * DS TWR RESP Main Loop (with STS enabled)
 *--------------------------------------------------------------------------- 
 */
void loop()
{
    // --- STS Key/IV Setup ---
    if (!firstLoopFlag)
    {
        // On first loop, load the STS key and IV into the chip.
        dwt_configurestskey(&cp_key);
        dwt_configurestsiv(&cp_iv);
        dwt_configurestsloadiv();
        firstLoopFlag = 1;
    }
    else
    {
        // On subsequent loops, reload the lower 32 bits of the STS IV (counter).
        dwt_writetodevice(STS_IV0_ID, 0, 4, (uint8_t *)&cp_iv);
        dwt_configurestsloadiv();
    }

    /* Reset preamble detection timeout and RX timeout for each new ranging process. */
    dwt_setpreambledetecttimeout(0);
    dwt_setrxtimeout(0);

    /* Activate reception immediately. */
    dwt_rxenable(DWT_START_RX_IMMEDIATE);

    /* Wait for a poll message or error/timeout. */
    while (!((status_reg = dwt_read32bitreg(SYS_STATUS_ID)) &
             (SYS_STATUS_RXFCG_BIT_MASK | SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR)))
    { /* spin */ };

    if (status_reg & SYS_STATUS_RXFCG_BIT_MASK)
    {
        /* Clear RX frame event flag. */
        dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_RXFCG_BIT_MASK);

        /* Read the received poll message into the buffer. */
        uint32_t frame_len = dwt_read32bitreg(RX_FINFO_ID) & FRAME_LEN_MAX_EX;
        if (frame_len <= RX_BUF_LEN)
        {
            dwt_readrxdata(rx_buffer, frame_len, 0);
        }

        /* Clear sequence number field for easier validation. */
        rx_buffer[ALL_MSG_SN_IDX] = 0;

        /* Validate that the received frame matches the expected poll message. */
        if (memcmp(rx_buffer, rx_poll_msg, ALL_MSG_COMMON_LEN) == 0)
        {
            uint32_t resp_tx_time;
            int ret;

            /* Retrieve poll reception timestamp. */
            poll_rx_ts = get_rx_timestamp_u64();

            /* Calculate the delayed transmission time for the response.
             * The delay is based on the poll reception timestamp plus a set delay.
             */
            resp_tx_time = (poll_rx_ts + (POLL_RX_TO_RESP_TX_DLY_UUS * UUS_TO_DWT_TIME)) >> 8;
            dwt_setdelayedtrxtime(resp_tx_time);

            /* Response TX timestamp is computed as the programmed time plus the TX antenna delay. */
            resp_tx_ts = (((uint64_t)(resp_tx_time & 0xFFFFFFFEUL)) << 8) + TX_ANT_DLY;

            /* Set RX delay and timeout for the final message reception. */
            dwt_setrxaftertxdelay(RESP_TX_TO_FINAL_RX_DLY_UUS);
            dwt_setrxtimeout(FINAL_RX_TIMEOUT_UUS);
            dwt_setpreambledetecttimeout(PRE_TIMEOUT);

            /* Prepare and send the response message. */
            tx_resp_msg[ALL_MSG_SN_IDX] = frame_seq_nb;
            dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_TXFRS_BIT_MASK);
            dwt_writetxdata(sizeof(tx_resp_msg), tx_resp_msg, 0);
            dwt_writetxfctrl(sizeof(tx_resp_msg), 0, 1);
            ret = dwt_starttx(DWT_START_TX_DELAYED | DWT_RESPONSE_EXPECTED);

            /* If starting the TX fails, skip this ranging exchange. */
            if (ret == DWT_ERROR)
            {
                return;
            }

            /* Wait for the final message (or timeout/error). */
            while (!((status_reg = dwt_read32bitreg(SYS_STATUS_ID)) &
                     (SYS_STATUS_RXFCG_BIT_MASK | SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR)))
            { /* spin */ };

            /* Increment the sequence number after sending the response. */
            frame_seq_nb++;

            if (status_reg & SYS_STATUS_RXFCG_BIT_MASK)
            {
                /* Clear RX and TX status flags. */
                dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_RXFCG_BIT_MASK | SYS_STATUS_TXFRS_BIT_MASK);

                /* Read the final message into the buffer. */
                frame_len = dwt_read32bitreg(RX_FINFO_ID) & FRAME_LEN_MAX_EX;
                if (frame_len <= RX_BUF_LEN)
                {
                    dwt_readrxdata(rx_buffer, frame_len, 0);
                }

                /* Clear sequence number field for validation. */
                rx_buffer[ALL_MSG_SN_IDX] = 0;
                if (memcmp(rx_buffer, rx_final_msg, ALL_MSG_COMMON_LEN) == 0)
                {
                    uint32_t poll_tx_ts, resp_rx_ts, final_tx_ts;
                    uint32_t poll_rx_ts_32, resp_tx_ts_32, final_rx_ts_32;
                    double Ra, Rb, Da, Db;
                    int64_t tof_dtu;

                    /* Retrieve response TX and final reception timestamps. */
                    resp_tx_ts = get_tx_timestamp_u64();
                    final_rx_ts = get_rx_timestamp_u64();

                    /* Extract the embedded timestamps from the final message. */
                    final_msg_get_ts(&rx_buffer[FINAL_MSG_POLL_TX_TS_IDX], &poll_tx_ts);
                    final_msg_get_ts(&rx_buffer[FINAL_MSG_RESP_RX_TS_IDX], &resp_rx_ts);
                    final_msg_get_ts(&rx_buffer[FINAL_MSG_FINAL_TX_TS_IDX], &final_tx_ts);

                    /* Compute the round-trip delays and then the time-of-flight. */
                    poll_rx_ts_32 = (uint32_t)poll_rx_ts;
                    resp_tx_ts_32 = (uint32_t)resp_tx_ts;
                    final_rx_ts_32 = (uint32_t)final_rx_ts;
                    Ra = (double)(resp_rx_ts - poll_tx_ts);
                    Rb = (double)(final_rx_ts_32 - resp_tx_ts_32);
                    Da = (double)(final_tx_ts - resp_rx_ts);
                    Db = (double)(resp_tx_ts_32 - poll_rx_ts_32);
                    tof_dtu = (int64_t)((Ra * Rb - Da * Db) / (Ra + Rb + Da + Db));

                    double tof = tof_dtu * DWT_TIME_UNITS;
                    double distance = tof * SPEED_OF_LIGHT;

                    // For example, print the distance (adjust as needed)
                    Serial.print(0); Serial.print(", ");
                    Serial.print(200); Serial.print(", ");
                    Serial.println(distance * 100);
                    
                    /* Optionally, wait a bit before starting the next exchange. */
                    Sleep(RNG_DELAY_MS);
                }
            }
            else
            {
                /* Clear RX error/timeout events if final message wasn't received. */
                dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR);
            }
        }
        else
        {
            /* Clear RX error/timeout events if poll message did not validate. */
            dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR);
        }
    }
    else
    {
        /* Clear RX error/timeout events if no poll message was received. */
        dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR);
    }
}
