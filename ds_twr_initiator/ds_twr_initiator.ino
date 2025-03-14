#include "dw3000.h"

// Application name
#define APP_NAME "DS TWR INIT v1.0 - DS with STS"

// Connection pins
const uint8_t PIN_RST = 15; // Reset pin
const uint8_t PIN_IRQ = 4;  // IRQ pin
const uint8_t PIN_SS  = 5;  // SPI select pin

/* Inter-ranging delay period, in milliseconds. */
#define RNG_DELAY_MS 500

/* Default antenna delay values for 64 MHz PRF. */
#define TX_ANT_DLY 16385
#define RX_ANT_DLY 16385

/* Frames used in the ranging process. */
static uint8_t tx_poll_msg[] = { 0x41, 0x88, 0, 0xCA, 0xDE, 'W', 'A', 'V', 'E', 0x21 };
static uint8_t rx_resp_msg[] = { 0x41, 0x88, 0, 0xCA, 0xDE, 'V', 'E', 'W', 'A', 0x10, 0x02, 0, 0 };
static uint8_t tx_final_msg[] = { 0x41, 0x88, 0, 0xCA, 0xDE, 'W', 'A', 'V', 'E', 0x23,
                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

/* Length of the common part of the message (first 10 bytes). */
#define ALL_MSG_COMMON_LEN 10

/* Indexes to access some fields in the frames. */
#define ALL_MSG_SN_IDX             2
#define FINAL_MSG_POLL_TX_TS_IDX  10
#define FINAL_MSG_RESP_RX_TS_IDX  14
#define FINAL_MSG_FINAL_TX_TS_IDX 18

/* Frame sequence number, incremented after each transmission. */
static uint8_t frame_seq_nb = 0;

/* Buffer to store the received response message. */
#define RX_BUF_LEN 20
static uint8_t rx_buffer[RX_BUF_LEN];

/* Status register copy for debugging. */
static uint32_t status_reg = 0;

/* Timing parameters (in UWB microseconds). */
#define POLL_TX_TO_RESP_RX_DLY_UUS 700
#define RESP_RX_TO_FINAL_TX_DLY_UUS 700
#define RESP_RX_TIMEOUT_UUS         300
#define PRE_TIMEOUT                 5

/* Timestamps for the ranging exchange. */
static uint64_t poll_tx_ts;
static uint64_t resp_rx_ts;
static uint64_t final_tx_ts;

/* External TX configuration (calibrated per device). */
extern dwt_txconfig_t txconfig_options;

/*
 * --- STS Key and IV Definitions ---
 * These 128-bit values (key and IV) must be identical on both the initiator and responder.
 * The defaults here are as specified in the IEEE 802.15.4z annex.
 */
static dwt_sts_cp_key_t cp_key = {
    0x14EB220F, 0xF86050A8, 0xD1D336AA, 0x14148674
};

static dwt_sts_cp_iv_t cp_iv = {
    0x1F9A3DE4, 0xD37EC3CA, 0xC44FA8FB, 0x362EEB34
};

/* Flag to ensure we load the STS key/IV once at startup. */
static uint8_t firstLoopFlag = 0;

/*
 * Modified configuration structure.
 * Note that the STS mode field is changed from DWT_STS_MODE_OFF to DWT_STS_MODE_AUTO
 * so that STS is enabled while preserving the DS three-message exchange.
 */
static dwt_config_t config = {
  5,                /* Channel number. */
  DWT_PLEN_128,     /* Preamble length (TX). */
  DWT_PAC8,         /* Preamble acquisition chunk size (RX). */
  9,                /* TX preamble code. */
  9,                /* RX preamble code. */
  1,                /* SFD type (non-standard 8-symbol SFD). */
  DWT_BR_6M8,       /* Data rate. */
  DWT_PHRMODE_STD,  /* PHY header mode. */
  DWT_PHRRATE_STD,  /* PHY header rate. */
  (129 + 8 - 8),    /* SFD timeout. */
  DWT_STS_MODE_1,/* Enable STS for DS TWR (was previously off). */
  DWT_STS_LEN_64,   /* STS length. */
  DWT_PDOA_M0       /* PDOA mode off. */
};

/* Function to send the poll message. */
static void send_tx_poll_msg(void)
{
    // Update sequence number in the poll message.
    tx_poll_msg[ALL_MSG_SN_IDX] = frame_seq_nb;
    
    // Clear previous TX frame sent status.
    dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_TXFRS_BIT_MASK);
    
    // Write poll message to the TX buffer.
    dwt_writetxdata(sizeof(tx_poll_msg), tx_poll_msg, 0);
    dwt_writetxfctrl(sizeof(tx_poll_msg), 0, 1);
    
    // Start immediate transmission.
    dwt_starttx(DWT_START_TX_IMMEDIATE);
    
    // Wait for TX frame sent event.
    while (!(dwt_read32bitreg(SYS_STATUS_ID) & SYS_STATUS_TXFRS_BIT_MASK))
    { /* spin */ }
    
    // Clear TX frame sent flag.
    dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_TXFRS_BIT_MASK);
}

void setup()
{
    UART_init();
    test_run_info((unsigned char *)APP_NAME);
    spiBegin(PIN_IRQ, PIN_RST);
    spiSelect(PIN_SS);
    
    /* Allow time for DW3000 startup (INIT_RC -> IDLE_RC). */
    Sleep(2);
    
    /* Ensure DW3000 is in IDLE state before proceeding. */
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
    
    /* Configure DW3000 with our STS-enabled DS configuration. */
    if (dwt_configure(&config))
    {
        UART_puts("CONFIG FAILED\r\n");
        while (1) ;
    }
    
    /* Configure TX spectrum parameters. */
    dwt_configuretxrf(&txconfig_options);
    
    /* Apply default antenna delay values. */
    dwt_setrxantennadelay(RX_ANT_DLY);
    dwt_settxantennadelay(TX_ANT_DLY);
    
    /* Set response delay, timeout, and preamble detection timeout. */
    dwt_setrxaftertxdelay(POLL_TX_TO_RESP_RX_DLY_UUS);
    dwt_setrxtimeout(RESP_RX_TIMEOUT_UUS);
    dwt_setpreambledetecttimeout(PRE_TIMEOUT);
    
    /* Enable LNA/PA and optionally LEDs for debugging. */
    dwt_setlnapamode(DWT_LNA_ENABLE | DWT_PA_ENABLE);
    dwt_setleds(DWT_LEDS_ENABLE | DWT_LEDS_INIT_BLINK);
}

void loop()
{
    // --- STS Key/IV Setup ---
    if (!firstLoopFlag)
    {
        // On first loop, configure and load the STS key and IV.
        dwt_configurestskey(&cp_key);
        dwt_configurestsiv(&cp_iv);
        dwt_configurestsloadiv();
        firstLoopFlag = 1;
    }
    else
    {
        // On subsequent loops, reload the lower 32 bits of the STS IV (the counter).
        dwt_writetodevice(STS_IV0_ID, 0, 4, (uint8_t *)&cp_iv);
        dwt_configurestsloadiv();
    }
    
    // --- DS TWR Exchange with STS Enabled ---
    
    // Send the poll message.
    tx_poll_msg[ALL_MSG_SN_IDX] = frame_seq_nb;
    dwt_writetxdata(sizeof(tx_poll_msg), tx_poll_msg, 0);
    dwt_writetxfctrl(sizeof(tx_poll_msg) + FCS_LEN, 0, 1);
    dwt_starttx(DWT_START_TX_IMMEDIATE | DWT_RESPONSE_EXPECTED);
    
    // Wait for response frame, timeout, or error.
    while (!( (status_reg = dwt_read32bitreg(SYS_STATUS_ID)) &
             (SYS_STATUS_RXFCG_BIT_MASK | SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR) ))
    { /* spin */ }
    
    // Increment frame sequence number.
    frame_seq_nb++;
    
    if (status_reg & SYS_STATUS_RXFCG_BIT_MASK)
    {
        // Clear RX and TX flags.
        dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_RXFCG_BIT_MASK | SYS_STATUS_TXFRS_BIT_MASK);
        
        // Read the received frame into the buffer.
        uint32_t frame_len = dwt_read32bitreg(RX_FINFO_ID) & FRAME_LEN_MAX_EX;
        if (frame_len <= RX_BUF_LEN)
        {
            dwt_readrxdata(rx_buffer, frame_len, 0);
        }
        
        // For validation, clear the sequence number field.
        rx_buffer[ALL_MSG_SN_IDX] = 0;
        
        // Validate that the response matches the expected header.
        if (memcmp(rx_buffer, rx_resp_msg, ALL_MSG_COMMON_LEN) == 0)
        {
            uint32_t final_tx_time;
            // Retrieve TX and RX timestamps.
            poll_tx_ts = get_tx_timestamp_u64();
            resp_rx_ts = get_rx_timestamp_u64();
            
            // Compute the final message transmission time.
            final_tx_time = (resp_rx_ts + (RESP_RX_TO_FINAL_TX_DLY_UUS * UUS_TO_DWT_TIME)) >> 8;
            dwt_setdelayedtrxtime(final_tx_time);
            
            // Compute the final TX timestamp (including antenna delay).
            final_tx_ts = (((uint64_t)(final_tx_time & 0xFFFFFFFEUL)) << 8) + TX_ANT_DLY;
            
            // Embed the timestamps in the final message.
            final_msg_set_ts(&tx_final_msg[FINAL_MSG_POLL_TX_TS_IDX], poll_tx_ts);
            final_msg_set_ts(&tx_final_msg[FINAL_MSG_RESP_RX_TS_IDX], resp_rx_ts);
            final_msg_set_ts(&tx_final_msg[FINAL_MSG_FINAL_TX_TS_IDX], final_tx_ts);
            
            // Send the final message.
            tx_final_msg[ALL_MSG_SN_IDX] = frame_seq_nb;
            dwt_writetxdata(sizeof(tx_final_msg), tx_final_msg, 0);
            dwt_writetxfctrl(sizeof(tx_final_msg) + FCS_LEN, 0, 1);
            
            int ret = dwt_starttx(DWT_START_TX_DELAYED);
            if (ret == DWT_SUCCESS)
            {
                while (!(dwt_read32bitreg(SYS_STATUS_ID) & SYS_STATUS_TXFRS_BIT_MASK))
                { /* spin */ }
                dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_TXFRS_BIT_MASK);
                frame_seq_nb++;
            }
        }
        else
        {
            // Response frame did not match expected header.
            UART_puts("BAD FRAME\r\n");
        }
    }
    else
    {
        // No good frame received; clear error/timeout flags.
        dwt_write32bitreg(SYS_STATUS_ID, SYS_STATUS_ALL_RX_TO | SYS_STATUS_ALL_RX_ERR | SYS_STATUS_TXFRS_BIT_MASK);
    }
    
    // Delay before the next ranging exchange.
    Sleep(RNG_DELAY_MS);
}
