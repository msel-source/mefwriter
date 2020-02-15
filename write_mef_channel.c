/**********************************************************************************************************************
 
 Copyright 2020, Mayo Foundation, Rochester MN. All rights reserved.
 
 This library contains functions to convert data samples to MEF version 3.0
 initialize_mef_channel_data() should be called first for each channel, which initializes the data in the channel
 structure.  Then write_mef_channel_data() is called with the actual sample data to be written to the mef.  Finally,
 close_mef_channel_file() will close out the channel mef file, and free allocated memory.
 
 To compile for a 64-bit intel system, linking with the following files is necessary:
 meflib.c, mefrec.c
 
 Usage and modification of this source code is governed by the Apache 2.0 license.
 You may not use this file except in compliance with this License.
 A copy of the Apache 2.0 License may be obtained at http://www.apache.org/licenses/LICENSE-2.0
 
 Thanks to all who acknowledge the Mayo Systems Electrophysiology Laboratory, Rochester, MN
 in academic publications of their work facilitated by this software.
 
 Written by Dan Crepeau 9/2011, originally to write MEF 2.1 files.
 
 Minor changes added 2/2012 to allow temporary index files to be written.
 Updated 10/2014 to write MEF 3.0.
 Updated 5/2015 for MEF 3.0 format changes.
 Updated 5/2016 for MEF 3.0 library updates and bug fixes.
 Updated 2/2017 to add support for exporting functions to .dll
 Updated 10/2017 initial support for MEF 3.0 annotations.
 
 When compiling for Windows .dll, define _EXPORT_FOR_DLL
 
 When using this code, the time-series data must be pre-sorted in increasing time order.  Unordered
 packets will cause data discontinuities.
 
 *************************************************************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
//#include <limits.h>

#include "write_mef_channel.h"

// TBD Won't be thread-safe without pthread, but for now ignore it
//pthread_mutex_t lock1 = PTHREAD_MUTEX_INITIALIZER;



// this creates a single static data block for use by all channels
// this assumes all channels will need the same size data block
ui1* GetDataBlockBuffer(ui8 block_len)
{
    static ui1 *buff = 0;
    
    if (buff == 0)
        buff = (ui1 *) malloc(block_len * 8);
    
    return buff;
}

#ifdef _EXPORT_FOR_DLL
__declspec(dllexport)
#endif

si4 append_mef_channel_data(CHANNEL_STATE *channel_state,
                            si1 *chan_map_name,
                            si4 new_segment_number,
                            si1 *mef_3_level_1_password,
                            si1 *mef_3_level_2_password,
                            si1 *mef3_session_directory,
                            ui8 num_secs_per_segment,
                            si4 bit_shift_flag
                            )
{
    extern MEF_GLOBALS	*MEF_globals;
    ui4 max_samps;
    SEGMENT     prev_segment;
    si1			prev_segment_name[MEF_SEGMENT_BASE_FILE_NAME_BYTES];
    si1         command[512];
    si1         extension[TYPE_BYTES];
    si1			mef3_session_path[MEF_FULL_FILE_NAME_BYTES], mef3_session_name[MEF_BASE_FILE_NAME_BYTES];
    si1			channel_path[MEF_FULL_FILE_NAME_BYTES], segment_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_SEGMENT_BASE_FILE_NAME_BYTES];
    UNIVERSAL_HEADER *uh;
    TIME_SERIES_METADATA_SECTION_2	*md2;
    METADATA_SECTION_3	*md3;
    
    char word[256];
    
    // must be a new segment greater than zero
    if (!(new_segment_number > 0))
        return 0;
    
    channel_state->if_appending = 1;
    
    sprintf(prev_segment_name, "%s/%s.timd/%s-%06d.segd", mef3_session_directory, chan_map_name, chan_map_name, new_segment_number - 1);
    
    (void) read_MEF_segment(&prev_segment, prev_segment_name, TIME_SERIES_CHANNEL_TYPE, mef_3_level_2_password, NULL, MEF_FALSE, MEF_FALSE);
    
    
    channel_state->raw_data_ptr_start = (si4 *) calloc((size_t) ((prev_segment.metadata_fps->metadata.time_series_section_2->block_interval / 1e6) *
                                                                 prev_segment.metadata_fps->metadata.time_series_section_2->sampling_frequency
                                                                 * 2), sizeof(si4));
    if (channel_state->raw_data_ptr_start == NULL)
    {
        fprintf(stderr, "Insufficient memory to allocate temporary channel buffer\n");
        exit(1);
    }
    channel_state->raw_data_ptr_current        = channel_state->raw_data_ptr_start;
    channel_state->block_hdr_time              = 0;
    channel_state->block_boundary              = 0;
    channel_state->last_chan_timestamp         = 0;
    channel_state->max_block_size              = 0;
    channel_state->max_block_len               = 0;
    channel_state->number_of_index_entries     = 0;
    channel_state->number_of_discontinuity_entries = 0;
    channel_state->block_sample_index          = 0;
    channel_state->number_of_samples           = 0;
    channel_state->discontinuity_flag          = 1;  // first block is by definition discontinuous
    channel_state->bit_shift_flag              = bit_shift_flag;
    channel_state->block_len                   = 0;  // this will be overwritten when write_mef_channel_data() is called  // TBD this variable isn't used
    
    channel_state->chan_num = prev_segment.metadata_fps->metadata.time_series_section_2->acquisition_channel_number;
    
    // get mef3 session name and path from passed directory
    extract_path_parts(mef3_session_directory, mef3_session_path, mef3_session_name, extension);
    MEF_snprintf(mef3_session_path, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", mef3_session_path, mef3_session_name, SESSION_DIRECTORY_TYPE_STRING);
    
    // make mef3 session directory
    // already exists
    //sprintf(command, "mkdir -p \"%s\" 2> /dev/null", mef3_session_path);
    //system(command);
    
    // set up a generic fps for universal header and password data
    channel_state->gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(channel_state->gen_fps, MEF_FALSE, MEF_FALSE, MEF_FALSE);
    uh = channel_state->gen_fps->universal_header;
    uh->segment_number = new_segment_number;
    MEF_strncpy(uh->session_name, prev_segment.metadata_fps->universal_header->session_name, MEF_BASE_FILE_NAME_BYTES);
    MEF_strncpy(uh->anonymized_name, prev_segment.metadata_fps->universal_header->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
    uh->start_time = UNIVERSAL_HEADER_START_TIME_NO_ENTRY;
    uh->end_time = UNIVERSAL_HEADER_END_TIME_NO_ENTRY;
    if (mef_3_level_2_password != NULL)
    {
        if (mef_3_level_1_password == NULL)
        {
            fprintf(stderr, "If a level 2 password is specified, then a level 1 password must be specified also.  Exiting...\n");
            exit(0);
        }
        channel_state->pwd = channel_state->gen_fps->password_data = process_password_data(NULL, mef_3_level_1_password, mef_3_level_2_password, uh);
    }
    else
        channel_state->pwd = channel_state->gen_fps->password_data = NULL;
    
    // make mef3 channel directory
    sprintf(channel_path, "%s/%s.%s", mef3_session_path, chan_map_name, TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING);
    // already exists
    //sprintf(command, "mkdir %s 2> /dev/null", channel_path);
    //system(command);
    
    // copy channel name into generic universal header
    MEF_strncpy(channel_state->gen_fps->universal_header->channel_name, prev_segment.metadata_fps->universal_header->channel_name, MEF_BASE_FILE_NAME_BYTES);
    
    // make mef3 segment name
    generate_segment_name(channel_state->gen_fps, segment_name);
    
    // save this for later for creating new segments
    strcpy(channel_state->channel_path, channel_path);
    
    // make mef3 segment directory
    sprintf(segment_path, "%s/%s.%s", channel_path, segment_name, SEGMENT_DIRECTORY_TYPE_STRING);
    sprintf(command, "mkdir %s", segment_path);
    system(command);
    
    // generate level UUID into generic universal_header
    generate_UUID(channel_state->gen_fps->universal_header->level_UUID);
    
    // set up mef3 time series metadata file
    channel_state->metadata_fps = allocate_file_processing_struct(METADATA_FILE_BYTES, TIME_SERIES_METADATA_FILE_TYPE_CODE, NULL, channel_state->gen_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(channel_state->metadata_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", segment_path, segment_name, TIME_SERIES_METADATA_FILE_TYPE_STRING);
    uh = channel_state->metadata_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = 1;
    uh->maximum_entry_size = METADATA_FILE_BYTES;
    initialize_metadata(channel_state->metadata_fps);
    if (mef_3_level_2_password != NULL)
    {
        channel_state->metadata_fps->metadata.section_1->section_2_encryption = LEVEL_1_ENCRYPTION_DECRYPTED;
        channel_state->metadata_fps->metadata.section_1->section_3_encryption = LEVEL_2_ENCRYPTION_DECRYPTED;
    }
    else
    {
        channel_state->metadata_fps->metadata.section_1->section_2_encryption = NO_ENCRYPTION;
        channel_state->metadata_fps->metadata.section_1->section_3_encryption = NO_ENCRYPTION;
    }
    md2 = channel_state->metadata_fps->metadata.time_series_section_2;
    MEF_strncpy(md2->channel_description, prev_segment.metadata_fps->metadata.time_series_section_2->channel_description, METADATA_CHANNEL_DESCRIPTION_BYTES);
    MEF_strncpy(md2->session_description, prev_segment.metadata_fps->metadata.time_series_section_2->session_description, METADATA_SESSION_DESCRIPTION_BYTES);
    md2->recording_duration = METADATA_RECORDING_DURATION_NO_ENTRY;
    md2->sampling_frequency = prev_segment.metadata_fps->metadata.time_series_section_2->sampling_frequency;
    md2->low_frequency_filter_setting = prev_segment.metadata_fps->metadata.time_series_section_2->low_frequency_filter_setting;
    md2->high_frequency_filter_setting = prev_segment.metadata_fps->metadata.time_series_section_2->high_frequency_filter_setting;
    md2->notch_filter_frequency_setting = prev_segment.metadata_fps->metadata.time_series_section_2->notch_filter_frequency_setting;
    md2->AC_line_frequency = prev_segment.metadata_fps->metadata.time_series_section_2->AC_line_frequency;
    md2->units_conversion_factor = prev_segment.metadata_fps->metadata.time_series_section_2->units_conversion_factor;
    MEF_strncpy(md2->units_description, "microvolts", TIME_SERIES_METADATA_UNITS_DESCRIPTION_BYTES);
    md2->maximum_native_sample_value = TIME_SERIES_METADATA_MAXIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY;  // must test against NaN later on
    md2->minimum_native_sample_value = TIME_SERIES_METADATA_MINIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY;
    md2->start_sample = prev_segment.metadata_fps->metadata.time_series_section_2->start_sample + prev_segment.metadata_fps->metadata.time_series_section_2->number_of_samples;
    md2->number_of_samples = 0;  // fill in when convert RED blocks
    md2->number_of_blocks = 0;  // fill in when convert RED blocks
    md2->maximum_block_bytes = 0;  // fill in when convert RED blocks
    md2->maximum_block_samples = 0;  // fill in when convert RED blocks
    md2->maximum_difference_bytes = 0;  // fill in when convert RED blocks
    md2->block_interval = prev_segment.metadata_fps->metadata.time_series_section_2->block_interval;
    md2->number_of_discontinuities = 0;  // fill in when convert RED blocks
    md2->maximum_contiguous_blocks = 0;  // fill in when convert RED blocks
    md2->maximum_contiguous_block_bytes = 0;  // fill in when convert RED blocks;
    md2->maximum_contiguous_samples = 0;  // fill in when convert RED blocks;
    md2->acquisition_channel_number = prev_segment.metadata_fps->metadata.time_series_section_2->acquisition_channel_number;  // for purposes of this program, these two will always be the same
    md3 = channel_state->metadata_fps->metadata.section_3;
    md3->recording_time_offset = prev_segment.metadata_fps->metadata.section_3->recording_time_offset;
    MEF_globals->recording_time_offset = md3->recording_time_offset;  // TBD make this thread_safe?
    md3->GMT_offset = prev_segment.metadata_fps->metadata.section_3->GMT_offset;
    MEF_globals->GMT_offset = md3->GMT_offset;  // TBD make this thread safe?
    //channel_state->gmt_offset_in_hours = gmt_offset;  // not used, since we already know offsets
    MEF_strncpy(md3->subject_name_1, prev_segment.metadata_fps->metadata.section_3->subject_name_1, METADATA_SUBJECT_NAME_BYTES);
    MEF_strncpy(md3->subject_name_2, prev_segment.metadata_fps->metadata.section_3->subject_name_2, METADATA_SUBJECT_NAME_BYTES);
    MEF_strncpy(md3->subject_ID, prev_segment.metadata_fps->metadata.section_3->subject_ID, METADATA_SUBJECT_ID_BYTES);
    MEF_strncpy(md3->recording_location, prev_segment.metadata_fps->metadata.section_3->recording_location, METADATA_RECORDING_LOCATION_BYTES);
    
    // set up mef3 time series indices file
    channel_state->ts_inds_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, TIME_SERIES_INDICES_FILE_TYPE_CODE, NULL, channel_state->metadata_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(channel_state->ts_inds_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", segment_path, segment_name, TIME_SERIES_INDICES_FILE_TYPE_STRING);
    uh = channel_state->ts_inds_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = 0;  // fill in when convert RED blocks
    uh->maximum_entry_size = TIME_SERIES_INDEX_BYTES;
    channel_state->ts_inds_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;  // write out the universal header, then the RED blocks piecemeal
    channel_state->ts_inds_fps->directives.close_file = MEF_FALSE;
    write_MEF_file(channel_state->ts_inds_fps);
    channel_state->inds_file_offset = UNIVERSAL_HEADER_BYTES;
    channel_state->ts_inds_fps->universal_header->body_CRC = CRC_START_VALUE;
    
    // set up mef3 time series data file
    channel_state->ts_data_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, TIME_SERIES_DATA_FILE_TYPE_CODE, NULL, channel_state->metadata_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(channel_state->ts_data_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", segment_path, segment_name, TIME_SERIES_DATA_FILE_TYPE_STRING);
    uh = channel_state->ts_data_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = 0;  // fill in when convert RED blocks
    uh->maximum_entry_size = 0;  // fill in when converet RED blocks
    channel_state->ts_data_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;  // write out the universal header, then the RED blocks piecemeal
    channel_state->ts_data_fps->directives.close_file = MEF_FALSE;
    write_MEF_file(channel_state->ts_data_fps);
    channel_state->data_file_offset = UNIVERSAL_HEADER_BYTES;
    channel_state->ts_data_fps->universal_header->body_CRC = CRC_START_VALUE;
    
    // allocate new memory for new RED blocks
    max_samps = (prev_segment.metadata_fps->metadata.time_series_section_2->block_interval / 1e6) *
    prev_segment.metadata_fps->metadata.time_series_section_2->sampling_frequency * 2;
    channel_state->rps = RED_allocate_processing_struct(max_samps, RED_MAX_COMPRESSED_BYTES(max_samps, 1), 0, RED_MAX_DIFFERENCE_BYTES(max_samps), 0, 0, channel_state->pwd);
    //channel_state->rps->directives.return_block_extrema = MEF_TRUE;
    // free original_data of rps, since we'll never use it anyway, and that pointer will be re-assigned before each RED compression
    // TBD modify it so original_data is not allocated in the first place?
    free(channel_state->rps->original_data);
    
    // set up discontinuity state information
    channel_state->discont_contiguous_blocks = 0;
    channel_state->discont_contiguous_samples = 0;
    channel_state->discont_contiguous_bytes = 0;
    
    // make these part of the channel state to keep everything thread-safe
    channel_state->out_data = (ui1 *) malloc(32000 * 8);  // This assumes 1 second blocks, sampled at 32000 Hz
    channel_state->temp_time_series_index = (ui1*) calloc(sizeof(ui1), 8+8+8+4+4+4+4+4+1+RED_BLOCK_PROTECTED_REGION_BYTES+RED_BLOCK_DISCRETIONARY_REGION_BYTES);
    
    channel_state->num_secs_per_segment = num_secs_per_segment;
    channel_state->next_segment_start_time = 0;
    channel_state->start_sample = 0;
    
    // TBD fix this
    //free_segment(&prev_segment, MEF_FALSE);
    
    return(1);
    
}

#ifdef _EXPORT_FOR_DLL
__declspec(dllexport)
#endif

si4 initialize_mef_channel_data ( CHANNEL_STATE *channel_state,
                                 sf8 secs_per_block,
                                 si1 *chan_map_name,
                                 si4 bit_shift_flag,
                                 sf8 low_frequency_filter_setting,
                                 sf8 high_frequency_filter_setting,
                                 sf8 notch_filter_frequency,
                                 sf8 AC_line_frequency,
                                 sf8 units_conversion_factor,
                                 si1 *channel_description,
                                 sf8 sampling_frequency,
                                 si8 block_interval,
                                 si4 chan_num,
                                 si1 *mef3_session_directory,
                                 sf4 gmt_offset,
                                 si1 *session_description,
                                 si1 *anonymized_subject_name,
                                 si1 *subject_first_name,
                                 si1 *subject_second_name,
                                 si1 *subject_id,
                                 si1 *institution,
                                 si1 *mef_3_level_1_password,
                                 si1 *mef_3_level_2_password,
                                 si1 *study_comments,
                                 si1 *channel_comments,
                                 ui8 num_secs_per_segment
                                 )
{
    extern int errno;
    extern MEF_GLOBALS	*MEF_globals;
    ui4 max_samps;
    si1	command[512];
    si1 extension[TYPE_BYTES];
    si1			mef3_session_path[MEF_FULL_FILE_NAME_BYTES], mef3_session_name[MEF_BASE_FILE_NAME_BYTES];
    si1			channel_path[MEF_FULL_FILE_NAME_BYTES], segment_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_SEGMENT_BASE_FILE_NAME_BYTES];
    UNIVERSAL_HEADER *uh;
    TIME_SERIES_METADATA_SECTION_2	*md2;
    METADATA_SECTION_3	*md3;
    
    channel_state->if_appending = 0;
    //fprintf(stderr, "sizeof struct: %d\n", sizeof(CHANNEL_STATE));
    
    //fprintf(stderr, "in initialize_mef_channel_data()\n");
    
    // add 10% to buffer size to account for possible sample frequency drift
    //fprintf(stderr,"%f, %f\n", secs_per_block, sampling_frequency);
    channel_state->raw_data_ptr_start = (si4 *) calloc((size_t) (secs_per_block * sampling_frequency * 2), sizeof(si4));
    if (channel_state->raw_data_ptr_start == NULL)
    {
        fprintf(stderr, "Insufficient memory to allocate temporary channel buffer\n");
        exit(1);
    }
    channel_state->raw_data_ptr_current        = channel_state->raw_data_ptr_start;
    channel_state->block_hdr_time              = 0;
    channel_state->block_boundary              = 0;
    channel_state->last_chan_timestamp         = 0;
    channel_state->max_block_size              = 0;
    channel_state->max_block_len               = 0;
    channel_state->number_of_index_entries     = 0;
    channel_state->number_of_discontinuity_entries = 0;
    channel_state->block_sample_index          = 0;
    channel_state->number_of_samples           = 0;
    channel_state->discontinuity_flag          = 1;  // first block is by definition discontinuous
    channel_state->bit_shift_flag              = bit_shift_flag;
    channel_state->block_len                   = 0;  // this will be overwritten when write_mef_channel_data() is called
    
    channel_state->chan_num = chan_num;
    
    
    // get mef3 session name and path from passed directory
    extract_path_parts(mef3_session_directory, mef3_session_path, mef3_session_name, extension);
    MEF_snprintf(mef3_session_path, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", mef3_session_path, mef3_session_name, SESSION_DIRECTORY_TYPE_STRING);
    
    //fprintf(stdout, "path: %s\n", mef3_session_path);
    // make mef3 session directory
    sprintf(command, "mkdir \"%s\" 2> /dev/null", mef3_session_path);
    system(command);
    
    // set up a generic fps for universal header and password data
    channel_state->gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(channel_state->gen_fps, MEF_FALSE, MEF_FALSE, MEF_FALSE);
    uh = channel_state->gen_fps->universal_header;
    uh->segment_number = 0;
    MEF_strncpy(uh->session_name, mef3_session_name, MEF_BASE_FILE_NAME_BYTES);
    MEF_strncpy(uh->anonymized_name, anonymized_subject_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
    uh->start_time = UNIVERSAL_HEADER_START_TIME_NO_ENTRY;
    uh->end_time = UNIVERSAL_HEADER_END_TIME_NO_ENTRY;
    if (mef_3_level_2_password != NULL)
    {
        if (mef_3_level_1_password == NULL)
        {
            fprintf(stderr, "If a level 2 password is specified, then a level 1 password must be specified also.  Exiting...\n");
            exit(0);
        }
        channel_state->pwd = channel_state->gen_fps->password_data = process_password_data(NULL, mef_3_level_1_password, mef_3_level_2_password, uh);
    }
    else
        channel_state->pwd = channel_state->gen_fps->password_data = NULL;
    
    // make mef3 channel directory
    sprintf(channel_path, "%s/%s.%s", mef3_session_path, chan_map_name, TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING);
    sprintf(command, "mkdir %s", channel_path);
    system(command);
    
    // copy channel name into generic universal header
    MEF_strncpy(channel_state->gen_fps->universal_header->channel_name, chan_map_name, MEF_BASE_FILE_NAME_BYTES);
    
    // make mef3 segment name
    generate_segment_name(channel_state->gen_fps, segment_name);
    
    // save this for later for creating new segments
    strcpy(channel_state->channel_path, channel_path);
    
    //fprintf(stdout, "segment path: %s\n", segment_path);
    // make mef3 segment directory
    sprintf(segment_path, "%s/%s.%s", channel_path, segment_name, SEGMENT_DIRECTORY_TYPE_STRING);
    sprintf(command, "mkdir %s", segment_path);
    system(command);
    
    // generate level UUID into generic universal_header
    generate_UUID(channel_state->gen_fps->universal_header->level_UUID);
    
    // set up mef3 time series metadata file
    channel_state->metadata_fps = allocate_file_processing_struct(METADATA_FILE_BYTES, TIME_SERIES_METADATA_FILE_TYPE_CODE, NULL, channel_state->gen_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(channel_state->metadata_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", segment_path, segment_name, TIME_SERIES_METADATA_FILE_TYPE_STRING);
    uh = channel_state->metadata_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = 1;
    uh->maximum_entry_size = METADATA_FILE_BYTES;
    initialize_metadata(channel_state->metadata_fps);
    if (mef_3_level_2_password != NULL)
    {
        channel_state->metadata_fps->metadata.section_1->section_2_encryption = LEVEL_1_ENCRYPTION_DECRYPTED;
        channel_state->metadata_fps->metadata.section_1->section_3_encryption = LEVEL_2_ENCRYPTION_DECRYPTED;
    }
    else
    {
        channel_state->metadata_fps->metadata.section_1->section_2_encryption = NO_ENCRYPTION;
        channel_state->metadata_fps->metadata.section_1->section_3_encryption = NO_ENCRYPTION;
    }
    md2 = channel_state->metadata_fps->metadata.time_series_section_2;
    MEF_strncpy(md2->channel_description, channel_comments, METADATA_CHANNEL_DESCRIPTION_BYTES);
    MEF_strncpy(md2->session_description, study_comments, METADATA_SESSION_DESCRIPTION_BYTES);
    md2->recording_duration = METADATA_RECORDING_DURATION_NO_ENTRY;
    md2->sampling_frequency = sampling_frequency;
    md2->low_frequency_filter_setting = low_frequency_filter_setting;
    md2->high_frequency_filter_setting = high_frequency_filter_setting;
    md2->notch_filter_frequency_setting = notch_filter_frequency;
    md2->AC_line_frequency = AC_line_frequency;
    md2->units_conversion_factor = units_conversion_factor;
    MEF_strncpy(md2->units_description, "microvolts", TIME_SERIES_METADATA_UNITS_DESCRIPTION_BYTES);
    md2->maximum_native_sample_value = TIME_SERIES_METADATA_MAXIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY;  // must test against NaN later on
    md2->minimum_native_sample_value = TIME_SERIES_METADATA_MINIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY;
    md2->start_sample = 0;
    md2->number_of_samples = 0;  // fill in when convert RED blocks
    md2->number_of_blocks = 0;  // fill in when convert RED blocks
    md2->maximum_block_bytes = 0;  // fill in when convert RED blocks
    md2->maximum_block_samples = 0;  // fill in when convert RED blocks
    md2->maximum_difference_bytes = 0;  // fill in when convert RED blocks
    md2->block_interval = block_interval;
    md2->number_of_discontinuities = 0;  // fill in when convert RED blocks
    md2->maximum_contiguous_blocks = 0;  // fill in when convert RED blocks
    md2->maximum_contiguous_block_bytes = 0;  // fill in when convert RED blocks;
    md2->maximum_contiguous_samples = 0;  // fill in when convert RED blocks;
    md2->acquisition_channel_number = chan_num;  // for purposes of this program, these two will always be the same
    md3 = channel_state->metadata_fps->metadata.section_3;
    md3->recording_time_offset = MEF_globals->recording_time_offset;
    md3->GMT_offset = MEF_globals->GMT_offset;  // TBD does this do anything?
    channel_state->gmt_offset_in_hours = gmt_offset;
    MEF_strncpy(md3->subject_name_1, subject_first_name, METADATA_SUBJECT_NAME_BYTES);
    MEF_strncpy(md3->subject_name_2, subject_second_name, METADATA_SUBJECT_NAME_BYTES);
    MEF_strncpy(md3->subject_ID, subject_id, METADATA_SUBJECT_ID_BYTES);
    MEF_strncpy(md3->recording_location, institution, METADATA_RECORDING_LOCATION_BYTES);
    
    // set up mef3 time series indices file
    channel_state->ts_inds_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, TIME_SERIES_INDICES_FILE_TYPE_CODE, NULL, channel_state->metadata_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(channel_state->ts_inds_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", segment_path, segment_name, TIME_SERIES_INDICES_FILE_TYPE_STRING);
    uh = channel_state->ts_inds_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = 0;  // fill in when convert RED blocks
    uh->maximum_entry_size = TIME_SERIES_INDEX_BYTES;
    channel_state->ts_inds_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;  // write out the universal header, then the RED blocks piecemeal
    channel_state->ts_inds_fps->directives.close_file = MEF_FALSE;
    write_MEF_file(channel_state->ts_inds_fps);
    channel_state->inds_file_offset = UNIVERSAL_HEADER_BYTES;
    channel_state->ts_inds_fps->universal_header->body_CRC = CRC_START_VALUE;
    
    // set up mef3 time series data file
    channel_state->ts_data_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, TIME_SERIES_DATA_FILE_TYPE_CODE, NULL, channel_state->metadata_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(channel_state->ts_data_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", segment_path, segment_name, TIME_SERIES_DATA_FILE_TYPE_STRING);
    uh = channel_state->ts_data_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = 0;  // fill in when convert RED blocks
    uh->maximum_entry_size = 0;  // fill in when converet RED blocks
    channel_state->ts_data_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;  // write out the universal header, then the RED blocks piecemeal
    channel_state->ts_data_fps->directives.close_file = MEF_FALSE;
    write_MEF_file(channel_state->ts_data_fps);
    channel_state->data_file_offset = UNIVERSAL_HEADER_BYTES;
    channel_state->ts_data_fps->universal_header->body_CRC = CRC_START_VALUE;
    
    // allocate new memory for new RED blocks
    max_samps = secs_per_block * sampling_frequency * 2;
    channel_state->rps = RED_allocate_processing_struct(max_samps, RED_MAX_COMPRESSED_BYTES(max_samps, 1), 0, RED_MAX_DIFFERENCE_BYTES(max_samps), 0, 0, channel_state->pwd);
    //channel_state->rps->directives.return_block_extrema = MEF_TRUE;
    // free original_data of rps, since we'll never use it anyway, and that pointer will be re-assigned before each RED compression
    // TBD modify it so original_data is not allocated in the first place?
    free(channel_state->rps->original_data);
    
    // set up discontinuity state information
    channel_state->discont_contiguous_blocks = 0;
    channel_state->discont_contiguous_samples = 0;
    channel_state->discont_contiguous_bytes = 0;
    
    // make these part of the channel state to keep everything thread-safe
    channel_state->out_data = (ui1 *) malloc(32000 * 8);  // This assumes 1 second blocks, sampled at 32000 Hz
    channel_state->temp_time_series_index = (ui1*) calloc(sizeof(ui1), 8+8+8+4+4+4+4+4+1+RED_BLOCK_PROTECTED_REGION_BYTES+RED_BLOCK_DISCRETIONARY_REGION_BYTES);
    
    channel_state->num_secs_per_segment = num_secs_per_segment;
    channel_state->next_segment_start_time = 0;
    channel_state->start_sample = 0;
    
    return(0);
}

si4 process_filled_block( CHANNEL_STATE *channel_state, si4* raw_data_ptr_start, ui4 num_entries,
                         ui8 block_len, si4 discontinuity_flag, ui8 block_hdr_time)
{
    ui1 *out_data;
    si4 *ddp, bit_shift_flag;
    ui8 i;
    extern MEF_GLOBALS	*MEF_globals;
    ui1 data_key[240];
    FILE *ofp;
    sf8 sampling_frequency;
    int chan_num;
    //static ui1 generated_offset = 0;
    //TIME_SERIES_INDEX time_series_index;
    ui1 *temp_time_series_index;
    sf8			temp_sf8;
    RED_PROCESSING_STRUCT	*rps;
    FILE_PROCESSING_STRUCT  *ts_inds_fps;
    FILE_PROCESSING_STRUCT  *ts_data_fps;
    FILE_PROCESSING_STRUCT  *metadata_fps;
    TIME_SERIES_METADATA_SECTION_2	*md2;
    UNIVERSAL_HEADER *uh_meta;
    UNIVERSAL_HEADER *uh_data;
    UNIVERSAL_HEADER *uh_inds;
    TIME_SERIES_INDEX temp_struct;
    
    
    // fprintf(stderr, "in process_filled_block\n");
    
    chan_num     = channel_state->chan_num;
    
    // bring in data from channel_state struct
    ofp                             = channel_state->out_file;
    bit_shift_flag                  = channel_state->bit_shift_flag;
    ts_data_fps                     = channel_state->ts_data_fps;
    ts_inds_fps                     = channel_state->ts_inds_fps;
    metadata_fps                    = channel_state->metadata_fps;
    temp_time_series_index          = channel_state->temp_time_series_index;
    rps                             = channel_state->rps;
    sampling_frequency              = channel_state->metadata_fps->metadata.time_series_section_2->sampling_frequency;
    
    // do nothing if there is nothing to be done
    if (num_entries == 0)
        return (0);
    if (block_len == 0)
        return (0);  // this should never happen, but check for it anyway
    
    //pthread_mutex_lock(&lock1);
    
    // only care about generating offset times if this is a brand-new session.
    // if we are appending to existing session, we alrady have offset times
    if ((channel_state->if_appending == 0) && (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT)))
    {
        // this only be done for one channel, assumes all channels have same offset
        if (MEF_globals->recording_time_offset == MEF_GLOBALS_RECORDING_TIME_OFFSET_DEFAULT) {
            // generate recording time offset & GMT offset
            //TBD bring in GMT offset
            generate_recording_time_offset(block_hdr_time, (si4) (channel_state->gmt_offset_in_hours * 3600.0));
            //generated_offset = 1;
        }
    }
    
    //pthread_mutex_unlock(&lock1);
    
    memset(data_key, 0, 240);  // for now, assume no data encryption
    
    // use a static out_data buffer, this buffer is shared across channels, so this will only
    // work if the block_len of all channels is the same.  If different channels need different
    // block_lens (ie, different sampling rates), then this block buffer will need to be part
    // of the channel struct.
    //out_data = GetDataBlockBuffer(block_len);
    
    // previous method (GetDataBlockBuffer) is not thread-safe, so use this method.
    // TBD clean this up, so it will works with any block_len, rather than hard-coding
    // 1 second blocks.
    out_data = channel_state->out_data;
    
    
    if (bit_shift_flag)
    {
        //shift 2 bits to 18 bit resolution
        ddp = raw_data_ptr_start;
        for(i = num_entries; i--;)
        {
            if (*ddp >= 0)
                *ddp++ = (si4) (((sf8) *ddp / (sf8) 4.0) + 0.5);
            else
                *ddp++ = (si4) (((sf8) *ddp / (sf8) 4.0) - 0.5);
        }
    }
    
    // set up RED compression
    rps->original_data = rps->original_ptr = raw_data_ptr_start;
    rps->directives.discontinuity = (discontinuity_flag== 1) ? MEF_TRUE : MEF_FALSE;
    rps->block_header->number_of_samples = num_entries;
    rps->block_header->start_time = block_hdr_time;
    
    // RED compress data block
    (void) RED_encode(rps);
    
    if (channel_state->num_secs_per_segment > 0 )
        check_for_new_segment(channel_state, rps->block_header->start_time);
    
    // write block to output file
    //fwrite(out_data, sizeof(si1), RED_block_size, ofp);
    (void) e_fwrite(rps->compressed_data, sizeof(ui1), (size_t) channel_state->rps->block_header->block_bytes, ts_data_fps->fp, ts_data_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    
    // update data body CRC
    ts_data_fps->universal_header->body_CRC = CRC_update(rps->compressed_data, rps->block_header->block_bytes, ts_data_fps->universal_header->body_CRC);
    
    // set recording_start_time on first pass
    uh_meta = channel_state->metadata_fps->universal_header;
    uh_data = channel_state->ts_data_fps->universal_header;
    uh_inds = channel_state->ts_inds_fps->universal_header;
    if (uh_meta->start_time == UNIVERSAL_HEADER_START_TIME_NO_ENTRY)
    {
        // needs to be offset, since universal header will always be written unencrypted
        // these timestamps are already offset, because the offsetting occurs in the RED compression routine
        uh_meta->start_time = rps->block_header->start_time;
        uh_data->start_time = rps->block_header->start_time;
        uh_inds->start_time = rps->block_header->start_time;
        
        // set start time for second segment (subtract because we're dealing with offset times)
        channel_state->next_segment_start_time = rps->block_header->start_time - (channel_state->num_secs_per_segment * 1e6);
    }
    
    md2 = metadata_fps->metadata.time_series_section_2;
    
    RED_find_extrema(rps->original_data, rps->block_header->number_of_samples, &temp_struct);
    
    // update segment  metadata files
    // maximum_native_sample_value
    if (md2->units_conversion_factor >= 0)
        temp_sf8 = (sf8) temp_struct.maximum_sample_value * md2->units_conversion_factor;
    else
        temp_sf8 = (sf8) temp_struct.minimum_sample_value * md2->units_conversion_factor;  // units_conversion_factor is negative, so use min value instead
    if (isnan(md2->maximum_native_sample_value))
        md2->maximum_native_sample_value = temp_sf8;
    if (temp_sf8 > md2->maximum_native_sample_value)
        md2->maximum_native_sample_value = temp_sf8;
    // minimum_native_sample_value
    if (md2->units_conversion_factor >= 0)
        temp_sf8 = (sf8) temp_struct.minimum_sample_value * md2->units_conversion_factor;
    else
        temp_sf8 = (sf8) temp_struct.maximum_sample_value * md2->units_conversion_factor;  // units_conversion_factor is negative, so use max value instead
    if (isnan(md2->minimum_native_sample_value))
        md2->minimum_native_sample_value = temp_sf8;
    if (temp_sf8 < md2->minimum_native_sample_value)
        md2->minimum_native_sample_value = temp_sf8;
    // maximum_block_bytes
    if (rps->block_header->block_bytes > md2->maximum_block_bytes)
        md2->maximum_block_bytes = rps->block_header->block_bytes;
    // maximum_difference_bytes
    if (rps->block_header->difference_bytes > md2->maximum_difference_bytes)
        md2->maximum_difference_bytes = rps->block_header->difference_bytes;
    // maximum_block_samples
    if (rps->block_header->number_of_samples > md2->maximum_block_samples)
        md2->maximum_block_samples = rps->block_header->number_of_samples;
    // number_of_samples
    md2->number_of_samples = md2->number_of_samples + num_entries;
    // number_of_blocks
    md2->number_of_blocks = md2->number_of_blocks + 1;
    // number_of_discontinuities
    if (discontinuity_flag == 1)
        md2->number_of_discontinuities = md2->number_of_discontinuities + 1;
    // in theory the next two only need to be set once, but we need to wait until this function, when we have real data,
    // in order to know what the offset and GMT times are.
    // recording_time_offset
    metadata_fps->metadata.section_3->recording_time_offset = MEF_globals->recording_time_offset;
    // GMT offset
    metadata_fps->metadata.section_3->GMT_offset = MEF_globals->GMT_offset;
    
    // update metadata recording_duration and end_time for all files
    uh_meta->end_time = block_hdr_time + (si8) (((((sf8) channel_state->rps->block_header->number_of_samples + 1) / md2->sampling_frequency) * (sf8) 1e6) + (sf8) 0.5);
    // needs to be offset, since universal header will always be written unencrypted
    if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
        apply_recording_time_offset(&uh_meta->end_time);
    uh_data->end_time = uh_meta->end_time;
    uh_inds->end_time = uh_meta->end_time;
    
    md2->recording_duration = uh_meta->end_time - uh_meta->start_time;
	// offset time values could be negative, so reverse sign if negative
	if (md2->recording_duration < 0)
		md2->recording_duration = 0 - md2->recording_duration;
    
    // update number_of_entries
    uh_data->number_of_entries = uh_data->number_of_entries + 1;
    uh_inds->number_of_entries = uh_inds->number_of_entries + 1;
    
    // update maximum_entry_size, this only applies to data, which is the largest number of samples in a RED block
    if (num_entries > uh_data->maximum_entry_size)
        uh_data->maximum_entry_size = num_entries;
    
    // set up block entry, TBD test to see if this works
    //time_series_index.file_offset = channel_state->data_file_offset;
    //time_series_index.start_time = channel_state->rps->block_header->start_time;
    //time_series_index.start_sample = start_sample;
    //time_series_index.number_of_samples = channel_state->rps->block_header->number_of_samples;
    //time_series_index.block_bytes = channel_state->rps->block_header->block_bytes;
    //time_series_index.maximum_sample_value = channel_state->rps->compression.maximum_sample_value;
    //time_series_index.minimum_sample_value = channel_state->rps->compression.minimum_sample_value;
    //time_series_index.flags =channel_state->rps->block_header->flags
    
    // write index entry to output mtf file, TBD, this might not be necessary
    memcpy(temp_time_series_index,    &(channel_state->data_file_offset),                          sizeof(ui8));
    memcpy(temp_time_series_index+8,  &(channel_state->rps->block_header->start_time),             sizeof(ui8));
    memcpy(temp_time_series_index+16, &(channel_state->start_sample),                                             sizeof(ui8));
    memcpy(temp_time_series_index+24, &(channel_state->rps->block_header->number_of_samples),      sizeof(ui4));
    memcpy(temp_time_series_index+28, &(channel_state->rps->block_header->block_bytes),            sizeof(ui4));
    memcpy(temp_time_series_index+32, &(temp_struct.maximum_sample_value),     sizeof(si4));
    memcpy(temp_time_series_index+36, &(temp_struct.minimum_sample_value),     sizeof(si4));
    memset(temp_time_series_index+40, 0, 4);
    memcpy(temp_time_series_index+44, &(channel_state->rps->block_header->flags),                  sizeof(ui1));
    
    // write block index entry
    //(void) e_fwrite(&block_index, sizeof(BLOCK_INDEX), (size_t) 1, block_inds_fps->fp, block_inds_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    (void) e_fwrite(temp_time_series_index, sizeof(ui1), (size_t) (45+RED_BLOCK_PROTECTED_REGION_BYTES+RED_BLOCK_DISCRETIONARY_REGION_BYTES), ts_inds_fps->fp, ts_inds_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    
    // update CRC
    ts_inds_fps->universal_header->body_CRC = CRC_update(temp_time_series_index, 45+RED_BLOCK_PROTECTED_REGION_BYTES+RED_BLOCK_DISCRETIONARY_REGION_BYTES, ts_inds_fps->universal_header->body_CRC);
    
    // update index file offset
    channel_state->inds_file_offset += (45+RED_BLOCK_PROTECTED_REGION_BYTES+RED_BLOCK_DISCRETIONARY_REGION_BYTES);
    
    // update discontinuity index
    if (discontinuity_flag == 1)
    {
        //channel_state->discont_block_number = number_of_index_entries;
        channel_state->discont_contiguous_blocks = 1;
        channel_state->discont_contiguous_samples = channel_state->rps->block_header->number_of_samples;
        channel_state->discont_contiguous_bytes = channel_state->rps->block_header->block_bytes;
        //channel_state->discont_contiguous_duration_start = block_hdr_time;
    }
    else
    {
        //channel_state->discont_block_number = doesn't change
        channel_state->discont_contiguous_blocks++;
        channel_state->discont_contiguous_samples += channel_state->rps->block_header->number_of_samples;
        channel_state->discont_contiguous_bytes += channel_state->rps->block_header->block_bytes;
    }
    
    // update metadata file
    // maximum_contiguous_blocks
    if (channel_state->discont_contiguous_blocks > md2->maximum_contiguous_blocks)
        md2->maximum_contiguous_blocks = channel_state->discont_contiguous_blocks;
    // maximum_contiguous_samples
    if (channel_state->discont_contiguous_samples > md2->maximum_contiguous_samples)
        md2->maximum_contiguous_samples = channel_state->discont_contiguous_samples;
    // maximum_contiguous_block_bytes
    if (channel_state->discont_contiguous_bytes > md2->maximum_contiguous_block_bytes)
        md2->maximum_contiguous_block_bytes = channel_state->discont_contiguous_bytes;
    
    // update fields for next time
    channel_state->data_file_offset += channel_state->rps->block_header->block_bytes;
    channel_state->start_sample += channel_state->rps->block_header->number_of_samples;
    
    // update mef header fields relating to block index
    channel_state->number_of_index_entries++;
    channel_state->number_of_samples += num_entries;
    
    // fprintf(stderr, "done with process_filled_block()");
    
    update_metadata(channel_state);  // necessary for real-time applications, otherwise, comment out this line.
    
    return(0);
}

#ifdef _EXPORT_FOR_DLL
__declspec (dllexport)
#endif

si4 initialize_meflib_dll()
{
    extern MEF_GLOBALS	*MEF_globals;
    
    (void)initialize_meflib();
 
    MEF_globals->recording_time_offset_mode = RTO_IGNORE;  // turn off timestamp offsetting by default
    
    return 0;
}

#ifdef _EXPORT_FOR_DLL
__declspec (dllexport)
#endif

si4 write_mef_channel_data( CHANNEL_STATE *channel_state,
                           ui8 *packet_times,
                           si4 *samps,
                           ui8 n_packets_to_process,
                           sf8 secs_per_block,
                           sf8 sampling_frequency)
{
    si4 *raw_data_ptr_start, *raw_data_ptr_current;
    ui8 block_len, block_hdr_time, block_boundary;
    ui8 last_chan_timestamp;
    si4 discontinuity_flag;
    si8 j;
    si8 block_interval;
    int chan_num;
    
    // fprintf(stderr, "in write_mef_channel_data()\n");
    
    chan_num          = channel_state->chan_num;
    
    // bring in data from channel_state struct
    raw_data_ptr_start   = channel_state->raw_data_ptr_start;
    raw_data_ptr_current = channel_state->raw_data_ptr_current;
    block_hdr_time       = channel_state->block_hdr_time;
    block_boundary       = channel_state->block_boundary;
    last_chan_timestamp  = channel_state->last_chan_timestamp;
    discontinuity_flag   = channel_state->discontinuity_flag;
    block_interval       = channel_state->metadata_fps->metadata.time_series_section_2->block_interval;
    
    // this is updated everytime, although it should never change between calls.
    // TBD add test to make sure it doesn't change?  This needs to be a parameter call because sometimes you don't
    // know the correct sampling frequency until data acutally arrives
    channel_state->metadata_fps->metadata.time_series_section_2->sampling_frequency = sampling_frequency;
    
    // set local constants
    block_len = (ui8) ceil(secs_per_block * sampling_frequency); //user-defined block size (s), convert to # of samples
    channel_state->block_len = block_len;
    
    for (j = 0; j < n_packets_to_process; ++j)
    {
        // set timestamp for the first block processed
        if (block_hdr_time == 0)
        {
            // block_hdr_time is the actual time put into the block header (timestamp of the first
            // block sample), while block_boundary is used only for calculation of which samples go
            // into which blocks.  block_boundary is never written to the mef file.
            block_hdr_time = packet_times[j];
            block_boundary = packet_times[j];
        }
        
        if ((abs((((si8)(packet_times[j]) - (si8)last_chan_timestamp))) >= DISCONTINUITY_TIME_THRESHOLD) ||
            (((si8)(packet_times[j]) - (si8)block_boundary) >= (si8)block_interval))
        {
            // Block needs to be compressed and written
            
            // See if data exists in the buffer before processing it.  Data might not exist if
            // this is the first sample we've processed so far.
            if ((raw_data_ptr_current - raw_data_ptr_start) > 0)
            {
                // process block of previously collected data
                process_filled_block(channel_state, raw_data_ptr_start, (raw_data_ptr_current - raw_data_ptr_start),
                                     block_len, discontinuity_flag, block_hdr_time);
            }
            
            // mark next block as being discontinuous if discontinuity is found
            if (abs((((si8)(packet_times[j]) - (si8)last_chan_timestamp))) >= DISCONTINUITY_TIME_THRESHOLD)
            {
                discontinuity_flag = 1;
                block_boundary = packet_times[j];
            }
            else
            {
                discontinuity_flag = 0;
                block_boundary += block_interval;
            }
            
            // set next block's timestamp
            block_hdr_time = packet_times[j];
            
            // move back to the beginning of the raw block
            raw_data_ptr_current = raw_data_ptr_start;
        }
        
        *raw_data_ptr_current++ = *samps++;
        
        last_chan_timestamp = packet_times[j];
    }
    
    // save state of channel for next time
    channel_state->raw_data_ptr_current = raw_data_ptr_current;
    channel_state->last_chan_timestamp  = last_chan_timestamp;
    channel_state->block_hdr_time       = block_hdr_time;
    channel_state->block_boundary       = block_boundary;
    channel_state->discontinuity_flag   = discontinuity_flag;
    
    return(0);
}

#ifdef _EXPORT_FOR_DLL
__declspec (dllexport)
#endif

si4 flush_mef_channel(CHANNEL_STATE *channel_state)
{
    si4 *raw_data_ptr_start, *raw_data_ptr_current;
    ui8 block_len, block_hdr_time, block_boundary;
    ui8 last_chan_timestamp;
    si4 discontinuity_flag;
    si8 block_interval;
    
    // bring in data from channel_state struct
    raw_data_ptr_start = channel_state->raw_data_ptr_start;
    raw_data_ptr_current = channel_state->raw_data_ptr_current;
    block_hdr_time = channel_state->block_hdr_time;
    block_boundary = channel_state->block_boundary;
    last_chan_timestamp = channel_state->last_chan_timestamp;
    discontinuity_flag = channel_state->discontinuity_flag;
    block_interval = channel_state->metadata_fps->metadata.time_series_section_2->block_interval;
    
    // this tests for the case where no data has yet been given to this channel.
    if (channel_state->block_len == 0)
        return 0;
    
    // See if data exists in the buffer before processing it.
    if ((raw_data_ptr_current - raw_data_ptr_start) > 0)
    {
        // process block of previously collected data
        process_filled_block(channel_state, raw_data_ptr_start, (raw_data_ptr_current - raw_data_ptr_start),
                             channel_state->block_len, discontinuity_flag, block_hdr_time);
    }
    
    // mark next block as being discontinuous
    discontinuity_flag = 1;
    
    // set to zero so it will be reset next packet
    block_boundary = 0;
    block_hdr_time = 0;
    
    // move back to the beginning of the raw block
    raw_data_ptr_current = raw_data_ptr_start;
    
    // save state of channel for next time
    channel_state->raw_data_ptr_current = raw_data_ptr_current;
    channel_state->last_chan_timestamp = last_chan_timestamp;
    channel_state->block_hdr_time = block_hdr_time;
    channel_state->block_boundary = block_boundary;
    channel_state->discontinuity_flag = discontinuity_flag;
    
    return (0);
}

si4 check_for_new_segment(CHANNEL_STATE *channel_state, ui8 start_time)
{
    
    FILE_PROCESSING_STRUCT  *ts_inds_fps;
    FILE_PROCESSING_STRUCT  *ts_data_fps;
    FILE_PROCESSING_STRUCT  *metadata_fps;
    UNIVERSAL_HEADER *uh;
    si1 segment_name[MEF_SEGMENT_BASE_FILE_NAME_BYTES];
    si1 segment_path[MEF_FULL_FILE_NAME_BYTES];
    si1	command[512];
    TIME_SERIES_METADATA_SECTION_2	*md2;
	extern MEF_GLOBALS	*MEF_globals;
    
    // ignore this function if we're still writing the first block to the first segment
    if (channel_state->next_segment_start_time == 0)
        return(0);
    
	if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
	{
		if (start_time > channel_state->next_segment_start_time)
			return(0);
	}
	else
	{
		if (start_time < channel_state->next_segment_start_time)
			return(0);
	}
    
    ts_inds_fps = channel_state->ts_inds_fps;
    ts_data_fps = channel_state->ts_data_fps;
    metadata_fps = channel_state->metadata_fps;
    
    // update and write segment metadata files as well as universal headers
    update_metadata(channel_state);
    
    // close old segment files
    fclose(ts_data_fps->fp);
    fclose(ts_inds_fps->fp);
    fclose(metadata_fps->fp);
    
    // set fp's to NULL, to force write_MEF_file() to do a new fopen()
    ts_data_fps->fp = NULL;
    ts_inds_fps->fp = NULL;
    metadata_fps->fp = NULL;
    
    // deal with data file
    uh = channel_state->ts_data_fps->universal_header;
    uh->segment_number++;
    generate_segment_name(ts_data_fps, segment_name);
    // make mef3 segment directory
    sprintf(segment_path, "%s/%s.%s", channel_state->channel_path, segment_name, SEGMENT_DIRECTORY_TYPE_STRING);
    sprintf(command, "mkdir %s", segment_path);
    //fprintf(stderr, "...%s...\n", command);
    system(command);
    // open new data file
    MEF_snprintf(ts_data_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", segment_path, segment_name, TIME_SERIES_DATA_FILE_TYPE_STRING);
    // update headers
    uh->start_time = start_time;
    uh->end_time = start_time;  // this will get overwritten very quickly
    uh->number_of_entries = 0;
    uh->maximum_entry_size = 0;
    generate_UUID(uh->level_UUID);
    generate_UUID(uh->file_UUID);
    ts_data_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;
    ts_data_fps->directives.close_file = MEF_FALSE;
    write_MEF_file(channel_state->ts_data_fps);
    channel_state->data_file_offset = UNIVERSAL_HEADER_BYTES;
    ts_data_fps->universal_header->body_CRC = CRC_START_VALUE;
    
    
    //deal with index file
    uh = channel_state->ts_inds_fps->universal_header;
    uh->segment_number++;
    // open new inds file
    MEF_snprintf(ts_inds_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", segment_path, segment_name, TIME_SERIES_INDICES_FILE_TYPE_STRING);
    // update headers
    uh->start_time = start_time;
    uh->end_time = start_time;  // this will get overwritten very quickly
    uh->number_of_entries = 0;
    uh->maximum_entry_size = TIME_SERIES_INDEX_BYTES;
    memcpy(uh->level_UUID, ts_data_fps->universal_header->level_UUID, 16);
    generate_UUID(uh->file_UUID);
    ts_inds_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;
    ts_inds_fps->directives.close_file = MEF_FALSE;
    write_MEF_file(channel_state->ts_inds_fps);
    channel_state->inds_file_offset = UNIVERSAL_HEADER_BYTES;
    ts_inds_fps->universal_header->body_CRC =  CRC_START_VALUE;
    
    //deal with metadata file
    uh = channel_state->metadata_fps->universal_header;
    uh->segment_number++;
    // open new inds file
    MEF_snprintf(metadata_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", segment_path, segment_name, TIME_SERIES_METADATA_FILE_TYPE_STRING);
    fps_open(metadata_fps, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    // update headers
    uh->body_CRC = CRC_START_VALUE;
    uh->start_time = start_time;
    uh->end_time = start_time;  // this will get overwritten very quickly
    uh->number_of_entries = 1;
    uh->maximum_entry_size = METADATA_FILE_BYTES;
    memcpy(uh->level_UUID, ts_data_fps->universal_header->level_UUID, 16);
    generate_UUID(uh->file_UUID);
    md2 = channel_state->metadata_fps->metadata.time_series_section_2;
    md2->recording_duration = METADATA_RECORDING_DURATION_NO_ENTRY;
    md2->maximum_native_sample_value = TIME_SERIES_METADATA_MAXIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY;  // must test against NaN later on
    md2->minimum_native_sample_value = TIME_SERIES_METADATA_MINIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY;
    md2->start_sample += md2->number_of_samples;  // start sample is incremented by previous segment's number_of_samples
    md2->number_of_samples = 0;
    md2->number_of_blocks = 0;
    md2->maximum_block_bytes = 0;
    md2->maximum_block_samples = 0;
    md2->maximum_difference_bytes = 0;
    md2->number_of_discontinuities = 0;
    md2->maximum_contiguous_blocks = 0;
    md2->maximum_contiguous_block_bytes = 0;
    md2->maximum_contiguous_samples = 0;
    
    // do internal channel_state variable resets
	if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
		channel_state->next_segment_start_time -= (channel_state->num_secs_per_segment * 1e6);
	else
		channel_state->next_segment_start_time += (channel_state->num_secs_per_segment * 1e6);
    // these internal variables are used to correctly update md2 values
    channel_state->discont_contiguous_blocks = 0;
    channel_state->discont_contiguous_samples = 0;
    channel_state->discont_contiguous_bytes = 0;
    // TBD are these fields still necessary?
    channel_state->number_of_index_entries = 0;
    channel_state->number_of_samples = 0;
    
    return(0);
}


#ifdef _EXPORT_FOR_DLL
// just a wrapper for update_metadata, since it needs to be called locally as
// well as exported to dll interface.
__declspec (dllexport) si4 update_metadata_dll(CHANNEL_STATE *channel_state)
{
    return update_metadata(channel_state);
}
#endif

si4 update_metadata(CHANNEL_STATE *channel_state)
{
    
    extern MEF_GLOBALS	*MEF_globals;
    int rewriting_metadata = 1;
    si1 *mode;
    struct stat	sb;
    
    // rewrite metadata file
    channel_state->metadata_fps->directives.close_file = MEF_FALSE;
    // this fseek might not be necessary, but it shouldn't hurt anything
    if (channel_state->metadata_fps->fp != NULL)
        e_fseek(channel_state->metadata_fps->fp, 0, SEEK_SET, channel_state->metadata_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    write_MEF_file(channel_state->metadata_fps);
    // decrypt if necessary, because write_MEF_file() encrypts if necessary
    // the "necessary" part is automatic in both cases
    decrypt_metadata(channel_state->metadata_fps);
    
    
    // re-calculate header CRC for index and data files.  Body CRCs for both files should already be up-to-date.
    channel_state->ts_inds_fps->universal_header->header_CRC = CRC_calculate(channel_state->ts_inds_fps->raw_data + CRC_BYTES, UNIVERSAL_HEADER_BYTES - CRC_BYTES);
    channel_state->ts_data_fps->universal_header->header_CRC = CRC_calculate(channel_state->ts_data_fps->raw_data + CRC_BYTES, UNIVERSAL_HEADER_BYTES - CRC_BYTES);
    
    // re-write data universal header and then go back to where we were
    e_fseek(channel_state->ts_data_fps->fp, 0, SEEK_SET, channel_state->ts_data_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    (void)e_fwrite(channel_state->ts_data_fps->universal_header, sizeof(UNIVERSAL_HEADER), (size_t)1, channel_state->ts_data_fps->fp, channel_state->ts_data_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    e_fseek(channel_state->ts_data_fps->fp, channel_state->data_file_offset, SEEK_SET, channel_state->ts_data_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    
    // re-write index universal header and then go back to where we were
    e_fseek(channel_state->ts_inds_fps->fp, 0, SEEK_SET, channel_state->ts_inds_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    (void)e_fwrite(channel_state->ts_inds_fps->universal_header, sizeof(UNIVERSAL_HEADER), (size_t)1, channel_state->ts_inds_fps->fp, channel_state->ts_inds_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    e_fseek(channel_state->ts_inds_fps->fp, channel_state->inds_file_offset, SEEK_SET, channel_state->ts_inds_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    
    // fprintf(stderr, "done update_metadata()\n");
    
    return(0);
}


#ifdef _EXPORT_FOR_DLL
__declspec (dllexport)
#endif

si4 close_mef_channel(CHANNEL_STATE *channel_state)
{
    
    // write remaining buffered data
    process_filled_block(channel_state,
                         channel_state->raw_data_ptr_start,
                         (channel_state->raw_data_ptr_current - channel_state->raw_data_ptr_start),
                         channel_state->block_len,
                         channel_state->discontinuity_flag,
                         channel_state->block_hdr_time);
    
    
    // update and write segment metadata files as well as universal headers
    update_metadata(channel_state);
    
    // close files
    fclose(channel_state->ts_data_fps->fp);
    fclose(channel_state->ts_inds_fps->fp);
    fclose(channel_state->metadata_fps->fp);
    
    return(0);
    
    // TBD skip freeing memory for now, c# isn't happy with it, may not even be necessary with garbage collection?
    
    // free memory
    free_file_processing_struct(channel_state->gen_fps);
    free_file_processing_struct(channel_state->ts_inds_fps);
    free_file_processing_struct(channel_state->metadata_fps);
    free_file_processing_struct(channel_state->ts_data_fps);
    
    free(channel_state->raw_data_ptr_start);
    channel_state->rps->original_data = NULL;  // original data was previously free'd and re-assigned, so no need to deallocate it here
    RED_free_processing_struct(channel_state->rps);
    free(channel_state->out_data);
    free(channel_state->temp_time_series_index);
    // TBD there appears to still be a small (368 byte) memory leak associated with each channel, it might be in meflib.c somewhere
    
    return(0);
}

#ifdef _EXPORT_FOR_DLL
__declspec (dllexport)
#endif

si4 create_or_append_annotations(ANNOTATION_STATE* annotation_state,
                                 si1* dir_name,
                                 sf4 gmt_offset,
                                 si1 *anonymized_subject_name)
{
    si1 file_name_temp[2048];
    FILE *file_temp;
    int file_exists;
    si1         extension[TYPE_BYTES];
    si1			mef3_session_path[MEF_FULL_FILE_NAME_BYTES], mef3_session_name[MEF_BASE_FILE_NAME_BYTES];
    
    //fprintf(stderr, "sizeof annotation_state struct: %d\n", sizeof(ANNOTATION_STATE));
    
    annotation_state->gmt_offset = gmt_offset;
    
    // set up a generic fps for universal header and password data
    annotation_state->gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(annotation_state->gen_fps, MEF_FALSE, MEF_FALSE, MEF_FALSE);
    
    file_exists = 0;
    
    extract_path_parts(dir_name, mef3_session_path, mef3_session_name, extension);
    
    // check to see if records file already exists
    sprintf(file_name_temp, "%s.mefd/%s.rdat", dir_name, mef3_session_name);
    if ((file_temp = fopen(file_name_temp, "r")))
    {
        fclose(file_temp);
        file_exists = 1;
    }
    
    if (file_exists)
    {
        FILE_PROCESSING_DIRECTIVES temp_directives;
        
        // set directives to read entire existing files, and keep file pointers open
        temp_directives.io_bytes = UNIVERSAL_HEADER_BYTES;
        temp_directives.close_file = MEF_FALSE;
        temp_directives.open_mode = FPS_R_OPEN_MODE;
        
        //fprintf(stderr, "Reading existing rdat and ridx files.\n");
        // read .rdat and .ridx files, leave file pointers open at the end of the files
        annotation_state->rdat_fps = read_MEF_file(NULL, file_name_temp, NULL, NULL, &(temp_directives), USE_GLOBAL_BEHAVIOR);
        sprintf(file_name_temp, "%s.mefd/%s.ridx", dir_name, mef3_session_name);
        annotation_state->ridx_fps = read_MEF_file(NULL, file_name_temp, NULL, NULL, &(temp_directives), USE_GLOBAL_BEHAVIOR);
        //fprintf(stderr, "Read existing rdat and ridx files.\n");
        
        // determine how big files are
        fseek(annotation_state->rdat_fps->fp, 0, SEEK_END);
        fseek(annotation_state->ridx_fps->fp, 0, SEEK_END);
        annotation_state->rdat_file_offset = (unsigned long) ftell(annotation_state->rdat_fps->fp);
        annotation_state->ridx_file_offset = (unsigned long) ftell(annotation_state->ridx_fps->fp);
        
        fclose(annotation_state->rdat_fps->fp);
        fclose(annotation_state->ridx_fps->fp);
    }
    else
    {
        // set up things common to both files' universal header
        MEF_strncpy(annotation_state->gen_fps->universal_header->session_name, mef3_session_name, MEF_BASE_FILE_NAME_BYTES);
        MEF_strncpy(annotation_state->gen_fps->universal_header->anonymized_name, anonymized_subject_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
        annotation_state->gen_fps->universal_header->start_time = UNIVERSAL_HEADER_START_TIME_NO_ENTRY;
        annotation_state->gen_fps->universal_header->end_time = UNIVERSAL_HEADER_END_TIME_NO_ENTRY;
        // generate level UUID into generic universal_header
        generate_UUID(annotation_state->gen_fps->universal_header->level_UUID);
        
        // allocate memory for new files
        annotation_state->rdat_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, RECORD_DATA_FILE_TYPE_CODE, NULL, annotation_state->gen_fps, UNIVERSAL_HEADER_BYTES);
        MEF_snprintf(annotation_state->rdat_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s.mefd/%s.rdat", dir_name, mef3_session_name);
        annotation_state->ridx_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, RECORD_INDICES_FILE_TYPE_CODE, NULL, annotation_state->gen_fps, UNIVERSAL_HEADER_BYTES);
        MEF_snprintf(annotation_state->ridx_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s.mefd/%s.ridx", dir_name, mef3_session_name);
        
        // create files (header only), leave file pointers open
        annotation_state->rdat_fps->directives.close_file = MEF_FALSE;
        annotation_state->rdat_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;
        annotation_state->rdat_fps->directives.open_mode = FPS_W_OPEN_MODE;
        annotation_state->ridx_fps->directives.close_file = MEF_FALSE;
        annotation_state->ridx_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;
        annotation_state->ridx_fps->directives.open_mode = FPS_W_OPEN_MODE;
        annotation_state->rdat_fps->universal_header->number_of_entries = 0;
        generate_UUID(annotation_state->rdat_fps->universal_header->file_UUID);
        write_MEF_file(annotation_state->rdat_fps);
        annotation_state->ridx_fps->universal_header->number_of_entries = 0;
        generate_UUID(annotation_state->ridx_fps->universal_header->file_UUID);
        write_MEF_file(annotation_state->ridx_fps);
        
        annotation_state->rdat_file_offset = UNIVERSAL_HEADER_BYTES;
        annotation_state->ridx_file_offset = UNIVERSAL_HEADER_BYTES;
        
        fclose(annotation_state->rdat_fps->fp);
        fclose(annotation_state->ridx_fps->fp);
    }
    
    return 0;
}

#ifdef _EXPORT_FOR_DLL
__declspec (dllexport)
#endif

si4 write_annotation(ANNOTATION_STATE* annotation_state,
                     ui8 unixTimestamp,
                     si1* type,
                     si4 code,
                     si1* annotation)
{
    extern MEF_GLOBALS	*MEF_globals;
    MEFREC_Seiz_1_0* mefrec_seiz;
    
    RECORD_HEADER *new_header;
    RECORD_INDEX *new_index;
    si4 pad_bytes;
    static const si1 pad_bytes_string[] = "~~~~~~~~~~~~~~~";  // 15 tildes, so we can fwrite between 0 and 15 of them to pad a record
    
    mefrec_seiz = NULL;
    
    if (!strcmp(type, "Siez") && !strcmp(type, "Note"))
        return 0;
    
    if (annotation_state->rdat_fps == NULL)
        return 0;
    
    if (annotation_state->ridx_fps == NULL)
        return 0;
    
    annotation_state->rdat_fps->fp = fopen(annotation_state->rdat_fps->full_file_name, "r+b");
    fseek(annotation_state->rdat_fps->fp, annotation_state->rdat_file_offset, SEEK_SET);
    annotation_state->ridx_fps->fp = fopen(annotation_state->ridx_fps->full_file_name, "r+b");
    fseek(annotation_state->ridx_fps->fp, annotation_state->ridx_file_offset, SEEK_SET);
    
    new_header = calloc(1, sizeof(RECORD_HEADER));
    new_index = calloc(1, sizeof(RECORD_INDEX));
    
    // populate header and index entry
    strcpy(new_header->type_string, type);
    strcpy(new_index->type_string, type);
    new_header->version_major = 1;
    new_index->version_major = 1;
    new_header->version_minor = 0;
    new_index->version_minor = 0;
    new_header->encryption = 0;
    new_index->encryption = 0;
    new_header->bytes = 0;
    
    // calculate bytes to put in header
    if (!strcmp(type, "Seiz"))
    {
        new_header->bytes = sizeof(MEFREC_Seiz_1_0);
    }
    if (!strcmp(type, "Note"))
    {
        new_header->bytes = strlen(annotation) + 1;  // add one for null terminator
    }
    
    // calculate pad bytes for possible encryption.  Encryption is done in 16 byte blocks.
    pad_bytes = 16 - (new_header->bytes % 16);
    if (pad_bytes == 16)
        pad_bytes = 0;
    
    new_header->bytes += pad_bytes;
    
	if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
	{
		// if we haven't already calculated an offset, do it now using our time and time zone
		if (MEF_globals->recording_time_offset == MEF_GLOBALS_RECORDING_TIME_OFFSET_DEFAULT) {
			generate_recording_time_offset(unixTimestamp, (si4)(annotation_state->gmt_offset * 3600.0));
		}
	}
    
    // these can be offset since they are not encrypted for both rdat and ridx
    new_header->time = unixTimestamp;
    if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
        apply_recording_time_offset(&(new_header->time));
    new_index->time = new_header->time;
    
    // update index offset before modifying rdat_file_offset variable, so we get beginning of record offset
    new_index->file_offset = annotation_state->rdat_file_offset;
    
    // keep track of where we are in rdat for next ridx entry
    annotation_state->rdat_file_offset += sizeof(RECORD_HEADER);
    
    long max_entry_size = 0;
    
    if (!strcmp(type, "Seiz"))
    {
        // populate seizure record body
        mefrec_seiz = calloc(1, sizeof(MEFREC_Seiz_1_0));
        sprintf(mefrec_seiz->annotation, annotation);
        mefrec_seiz->earliest_onset = unixTimestamp;
        mefrec_seiz->latest_offset = unixTimestamp;
        mefrec_seiz->duration = 0;
        mefrec_seiz->number_of_channels = 4;
        mefrec_seiz->onset_code = code;
        
        // keep track of where we are in rdat for next ridx entry
        annotation_state->rdat_file_offset += sizeof(MEFREC_Seiz_1_0);
        
        // keep track of entry size to update universal header fro rdat
        max_entry_size = sizeof(MEFREC_Seiz_1_0) + sizeof(RECORD_HEADER);
        
        // calculate CRC
        new_header->record_CRC = CRC_calculate((ui1*)new_header + CRC_BYTES, sizeof(RECORD_HEADER) - CRC_BYTES);
        new_header->record_CRC = CRC_update((ui1*)mefrec_seiz, (size_t)sizeof(MEFREC_Seiz_1_0), new_header->record_CRC);
    }
    if (!strcmp(type, "Note"))
    {
        // calculate offset for record data file
        annotation_state->rdat_file_offset += strlen(annotation) + 1;  // add one for null terminator
        
        // calculate offset for record data file
        max_entry_size = strlen(annotation) + 1 + sizeof(RECORD_HEADER);
        
        // calculate CRC
        new_header->record_CRC = CRC_calculate((ui1*)new_header + CRC_BYTES, sizeof(RECORD_HEADER) - CRC_BYTES);
        new_header->record_CRC = CRC_update((ui1*)annotation, (size_t)(strlen(annotation) + 1), new_header->record_CRC);
    }
    
    // account for pad bytes, this code will be the same regardless of record type
    annotation_state->rdat_file_offset += pad_bytes;
    max_entry_size += pad_bytes;
    if (pad_bytes)
        new_header->record_CRC = CRC_update((ui1*)pad_bytes_string, (size_t)pad_bytes, new_header->record_CRC);
    
    // we know the CRC for the record header (crc of header and body) so now we can write them.
    (void)e_fwrite(new_header, sizeof(ui1), (size_t)sizeof(RECORD_HEADER), annotation_state->rdat_fps->fp, annotation_state->rdat_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    if (!strcmp(type, "Seiz"))
    {
        (void)e_fwrite(mefrec_seiz, sizeof(ui1), (size_t)sizeof(MEFREC_Seiz_1_0), annotation_state->rdat_fps->fp, annotation_state->rdat_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    }
    if (!strcmp(type, "Note"))
    {
        (void)e_fwrite(annotation, sizeof(ui1), (size_t)(strlen(annotation) + 1), annotation_state->rdat_fps->fp, annotation_state->rdat_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    }
    // write pad bytes, if necessary
    if (pad_bytes)
    {
        (void)e_fwrite(pad_bytes_string, sizeof(si1), (size_t)pad_bytes, annotation_state->rdat_fps->fp, annotation_state->rdat_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    }
    // write index
    (void)e_fwrite(new_index, sizeof(ui1), (size_t)sizeof(RECORD_INDEX), annotation_state->ridx_fps->fp, annotation_state->ridx_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    annotation_state->ridx_file_offset += sizeof(RECORD_INDEX);
    
    /*
     // done appending, close files
     fclose(annotation_state->rdat_fps->fp);
     fclose(annotation_state->ridx_fps->fp);
     
     // reopen files, so we can modify the header
     annotation_state->rdat_fps->fp = fopen(annotation_state->rdat_fps->full_file_name, "w+b");
     fseek(annotation_state->rdat_fps->fp, 0, SEEK_SET);
     annotation_state->ridx_fps->fp = fopen(annotation_state->ridx_fps->full_file_name, "w+b");
     fseek(annotation_state->ridx_fps->fp, 0, SEEK_SET);
     */
    
    // udpdate universal_header fields
    // update start_time, if necessary
    if (annotation_state->rdat_fps->universal_header->start_time == UNIVERSAL_HEADER_START_TIME_NO_ENTRY)
    {
        annotation_state->rdat_fps->universal_header->start_time = unixTimestamp;
        annotation_state->ridx_fps->universal_header->start_time = unixTimestamp;
        // apply offset, since universal header is always written unencrypted
        if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
        {
            apply_recording_time_offset(&annotation_state->rdat_fps->universal_header->start_time);
            apply_recording_time_offset(&annotation_state->ridx_fps->universal_header->start_time);
        }
    }
    
    // update end_time
    annotation_state->rdat_fps->universal_header->end_time = unixTimestamp;
    annotation_state->ridx_fps->universal_header->end_time = unixTimestamp;
    // apply offset, since universal header is always written unencrypted
    if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
    {
        apply_recording_time_offset(&annotation_state->rdat_fps->universal_header->end_time);
        apply_recording_time_offset(&annotation_state->ridx_fps->universal_header->end_time);
    }
    
    // update max_entry_size, if necessary
    if ((annotation_state->rdat_fps->universal_header->maximum_entry_size < max_entry_size) ||
        (annotation_state->rdat_fps->universal_header->maximum_entry_size == UNIVERSAL_HEADER_MAXIMUM_ENTRY_SIZE_NO_ENTRY))
    {
        annotation_state->rdat_fps->universal_header->maximum_entry_size = max_entry_size;
        annotation_state->ridx_fps->universal_header->maximum_entry_size = max_entry_size;
    }
    
    
    // update body CRCs
    annotation_state->rdat_fps->universal_header->body_CRC = CRC_update((ui1*)new_header, (size_t)sizeof(RECORD_HEADER), annotation_state->rdat_fps->universal_header->body_CRC);
    if (!strcmp(type, "Seiz"))
    {
        annotation_state->rdat_fps->universal_header->body_CRC = CRC_update((ui1*)mefrec_seiz, (size_t)sizeof(MEFREC_Seiz_1_0), annotation_state->rdat_fps->universal_header->body_CRC);
    }
    annotation_state->ridx_fps->universal_header->body_CRC = CRC_update((ui1*)new_index,  (size_t)sizeof(RECORD_INDEX),  annotation_state->ridx_fps->universal_header->body_CRC);
    
    // update number_of_entries for both files
    annotation_state->rdat_fps->universal_header->number_of_entries++;
    annotation_state->ridx_fps->universal_header->number_of_entries++;
    
    // re-calculate header CRC for index and data files.  Body CRCs for both files should already be up-to-date.
    annotation_state->rdat_fps->universal_header->header_CRC = CRC_calculate((ui1*)annotation_state->rdat_fps->raw_data + CRC_BYTES, UNIVERSAL_HEADER_BYTES - CRC_BYTES);
    annotation_state->ridx_fps->universal_header->header_CRC = CRC_calculate((ui1*)annotation_state->ridx_fps->raw_data + CRC_BYTES, UNIVERSAL_HEADER_BYTES - CRC_BYTES);
    
    // rewrite universal headers
    // re-write data universal header and then go back to where we were
    e_fseek(annotation_state->rdat_fps->fp, 0, SEEK_SET, annotation_state->rdat_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    (void)e_fwrite(annotation_state->rdat_fps->universal_header, sizeof(UNIVERSAL_HEADER), (size_t)1, annotation_state->rdat_fps->fp, annotation_state->rdat_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    e_fseek(annotation_state->rdat_fps->fp, annotation_state->rdat_file_offset, SEEK_SET, annotation_state->rdat_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    
    // re-write data universal header and then go back to where we were
    e_fseek(annotation_state->ridx_fps->fp, 0, SEEK_SET, annotation_state->ridx_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    (void)e_fwrite(annotation_state->ridx_fps->universal_header, sizeof(UNIVERSAL_HEADER), (size_t)1, annotation_state->ridx_fps->fp, annotation_state->ridx_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    e_fseek(annotation_state->ridx_fps->fp, annotation_state->ridx_file_offset, SEEK_SET, annotation_state->ridx_fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    
    // free memory
    if (!strcmp(type, "Seiz"))
    {
        free(mefrec_seiz);
    }
    free(new_header);
    free(new_index);
    
    return 0;
}

#ifdef _EXPORT_FOR_DLL
__declspec (dllexport)
#endif 

si4 close_annotation(ANNOTATION_STATE* annotation_state)
{
    // TBD clean up memory
    
    // close files
    if (annotation_state->rdat_fps->fp != NULL)
        fclose(annotation_state->rdat_fps->fp);
    if (annotation_state->ridx_fps->fp != NULL)
        fclose(annotation_state->ridx_fps->fp);
    
    return 0;
}
