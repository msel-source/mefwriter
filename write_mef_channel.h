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


#ifndef WRITE_MEF_CHANNEL
#define WRITE_MEF_CHANNEL

#ifdef __cplusplus
extern "C" {
#endif
    
#include "meflib.h"
    
    typedef struct SESSION_STATE SESSION_STATE;
    
    typedef struct {
        si4     chan_num;
        RED_PROCESSING_STRUCT	*rps;
        si4     *raw_data_ptr_start;
        si4     *raw_data_ptr_current;
        ui8     block_hdr_time;
        ui8     block_boundary;
        ui8     last_chan_timestamp;
        ui8     max_block_size;
        ui8     max_block_len;
        ui8     number_of_index_entries;
        ui8     number_of_discontinuity_entries;
        ui8     number_of_samples;
        ui8     block_sample_index;
        si8     start_sample;
        si4     discontinuity_flag;
        si4     bit_shift_flag;
        FILE    *out_file;
        ui1*    out_data;
        ui1*    temp_time_series_index;
        ui4     discont_contiguous_blocks;
        si8     discont_contiguous_samples;
        si8     discont_contiguous_bytes;
        ui8     block_len;
        FILE_PROCESSING_STRUCT *gen_fps;
        FILE_PROCESSING_STRUCT *metadata_fps;
        FILE_PROCESSING_STRUCT *ts_data_fps;
        FILE_PROCESSING_STRUCT *ts_inds_fps;
        si8     inds_file_offset;
        si8     data_file_offset;
        PASSWORD_DATA		*pwd;
        sf4 gmt_offset_in_hours;
        si1 segment_name[MEF_SEGMENT_BASE_FILE_NAME_BYTES];
        si1 segment_path[MEF_FULL_FILE_NAME_BYTES];
        ui8 num_secs_per_segment;
        ui8 next_segment_start_time;
        si1 channel_path[MEF_FULL_FILE_NAME_BYTES];
        si1    if_appending;
    } CHANNEL_STATE;
    
    typedef struct {
        FILE_PROCESSING_STRUCT *gen_fps;
        FILE_PROCESSING_STRUCT *rdat_fps;
        FILE_PROCESSING_STRUCT *ridx_fps;
        sf4 gmt_offset;
        si8     rdat_file_offset;
        si8     ridx_file_offset;
    } ANNOTATION_STATE;
    
    // Subroutine declarations
    
#ifndef _EXPORT_FOR_DLL
     si4 initialize_mef_channel_data(CHANNEL_STATE *channel_state,
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
     si1 *subject_ID,
     si1 *institution,
     si1 *mef_3_level_1_password,
     si1 *mef_3_level_2_password,
     si1 *study_comments,
     si1 *channel_comments,
     ui8 num_secs_per_segment
     );
#endif
    
#ifndef _EXPORT_FOR_DLL
     si4 write_mef_channel_data(CHANNEL_STATE *channel_state,
     ui8 *packet_times,
     si4 *samps,
     ui8 n_packets_to_process,
     sf8 secs_per_block,
     sf8 sampling_frequency);
#endif
    si4 check_for_new_segment(CHANNEL_STATE *channel_state, ui8 start_time);
    
    si4 update_metadata(CHANNEL_STATE *channel_state);

#ifndef _EXPORT_FOR_DLL
     si4 close_mef_channel(CHANNEL_STATE *channel_state);
     
     si4 append_mef_channel_data(CHANNEL_STATE *channel_state,
     si1 *chan_map_name,
     si4 new_segment_number,
     si1 *mef_3_level_1_password,
     si1 *mef_3_level_2_password,
     si1 *mef3_session_directory,
     ui8 num_secs_per_segment,
     si4 bit_shift_flag
     );
#endif
    
#ifndef _EXPORT_FOR_DLL
    si4 create_or_append_annotations(ANNOTATION_STATE* annotation_state,
                                     si1* dir_name,
                                     sf4 gmt_offset,
                                     si1 *anonymized_subject_name);
    si4 write_annotation(ANNOTATION_STATE* annotation_state,
                         ui8 unixTimestamp,
                         si1* type,
                         void* record);
    si4 close_annotation(ANNOTATION_STATE* annotation_state);
#endif
    
    

#define DISCONTINUITY_TIME_THRESHOLD 100000   // 100000 microseconds = .1 seconds
    
    
#ifdef __cplusplus
}
#endif 

#endif
