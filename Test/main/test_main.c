#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "../../Include/dr_wav.h"
#include "../../Include/getopt.h"
#include "../../Include/logger.h"
#include "../../Include/ini.h"

#define FRAME_SIZE                      512
#define FRAME_MOVE                      256
#define MAX_CHANNEL                     16
#define MAX_CHANNEL_SAMPLE              FRAME_MOVE * MAX_CHANNEL
#define PATH_LEN                        1024
#define FS                              16000
#define MIC_NUM                         2

short in_audio[MAX_CHANNEL_SAMPLE];
short out_audio[MIC_NUM * FRAME_MOVE];
drwav in_wav, out_wav;
char config_filename[PATH_LEN] = { 0 }, 
	 in_wav_filename[PATH_LEN] = { 0 },
     out_wav_filename[PATH_LEN] = { 0 },
	 log_filename[PATH_LEN] = {0};

void parse_command_line(int argc, char* argv[])
{
	int oc = 0;
	while ((oc = getopt(argc, argv, "i:o:c:l:h")) != -1) {
		switch (oc) {
		case 'i':
			strcpy(in_wav_filename, optarg);
			break;
		case 'o':
			strcpy(out_wav_filename, optarg);
			break;
		case 'c':
			strcpy(config_filename, optarg);
			break;
		case 'l':
			strcpy(log_filename, optarg);
			break;
		case 'h':
			return;
		default:
			break;
		}
	}
}

typedef struct
{
	int version;
	const char name[PATH_LEN];
	const char email[PATH_LEN];
} configuration;

static int parse_handler(void* user, const char* section, const char* name,
	const char* value)
{
	configuration* pconfig = (configuration*)user;
#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
	if (MATCH("protocol", "version")) {
		pconfig->version = atoi(value);
	}
	else if (MATCH("user", "name")) {
		strcpy(pconfig->name, value);
	}
	else if (MATCH("user", "email")) {
		strcpy(pconfig->email, value);
	}
	else {
		return 0;  /* unknown section/name, error */
	}

	return 1;
}

void main(int argc, char* argv[])
{
	char* argk[] = { " ",
		"-i","./data/test.wav",
		"-o","./data/test_out.wav",
		"-l","./data/test.log",
		"-c","./data/test.ini"};

	if (argc < 4) {
		argc = sizeof(argk) / sizeof(argk[0]);
		argv = argk;
	}

	parse_command_line(argc, argv);
	logger_initFileLogger(log_filename, 1024 * 1024, 5);
	logger_setLevel(LogLevel_DEBUG);
	LOG_INFO("input file name:%s", in_wav_filename);
	LOG_INFO("output file name:%s", out_wav_filename);
	LOG_INFO("config file name:%s", config_filename);
	LOG_INFO("log file name:%s", log_filename);

	configuration config;
	if (ini_parse(config_filename, parse_handler, &config) < 0) {
		LOG_ERROR("Can't load %s", config_filename);
		return;
	}

	if (!drwav_init_file(&in_wav, in_wav_filename, NULL)) {
		LOG_ERROR("Error opening WAV file:%s", in_wav_filename);
		return;
	}

	long flen = in_wav.totalPCMFrameCount;
	long n_samples = 0;
	drwav_data_format format;
	format.container = drwav_container_riff;     // <-- drwav_container_riff = normal WAV files, drwav_container_w64 = Sony Wave64.
	format.format = DR_WAVE_FORMAT_PCM;          // <-- Any of the DR_WAVE_FORMAT_* codes.
	format.channels = 1;
	format.sampleRate = FS;
	format.bitsPerSample = 16;
	drwav_init_file_write(&out_wav, out_wav_filename, &format, NULL);
	
	float in_data[FRAME_SIZE + 2], out_data[FRAME_SIZE];
	
	while (flen > 0) {
		n_samples = drwav_read_pcm_frames_s16(&in_wav, FRAME_MOVE, in_audio);

		drwav_uint64 framesWritten = drwav_write_pcm_frames(&out_wav, n_samples, in_audio);

		flen -= n_samples;
	}

	drwav_uninit(&in_wav);
	drwav_uninit(&out_wav);
}