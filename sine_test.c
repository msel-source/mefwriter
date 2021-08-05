// Multiscale Electrophysiology Format (MEF) version 3.0
// Copyright 2021, Mayo Foundation, Rochester MN. All rights reserved.
// Written by Matt Stead, Ben Brinkmann, and Dan Crepeau.

// Usage and modification of this source code is governed by the Apache 2.0 license.
// You may not use this file except in compliance with this License.
// A copy of the Apache 2.0 License may be obtained at http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under this License is distributed on an "as is" basis,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Thanks to all who acknowledge the Mayo Systems Electrophysiology Laboratory, Rochester, MN
// in academic publications of their work facilitated by this software.

// This code is a simple example in C of how to create a MEF 3.0 channel, add some data to it,
// then close the channel.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "write_mef_channel.h"

extern MEF_GLOBALS	*MEF_globals;

int main()
{
    int i;
    sf8 sampling_frequency;
    sf8 sine_amplitude;
    sf8 sine_frequency;
    sf8 seconds_per_block;
    CHANNEL_STATE *mef_channel_state_struct;
    ANNOTATION_STATE *annotation_state_struct;
    si4 *samps;
    ui8 *packet_times;
    si8 base_timestamp;
    char dir_name[128];
    char temp_name[MEF_FULL_FILE_NAME_BYTES];
    MEFREC_Curs_1_0 *cursor_pointer;
    MEFREC_Epoc_1_0 *epoch_pointer;
    FILE_PROCESSING_STRUCT *records_fps;
    
    // initialize MEF3 library
    (void) initialize_meflib();
    MEF_globals->recording_time_offset_mode = RTO_IGNORE;  // turn off timestamp offsetting by default
    
    // initialize variables
    sampling_frequency = 1000.0;   // Hz
    seconds_per_block = 1.0;
#ifdef _WIN32
    sprintf(dir_name, "c:\\sine_test");
#else
    sprintf(dir_name, "sine_test");
#endif
    
    // allocate buffers
    samps = (si4*) calloc((size_t) 10000, sizeof(si4));
    packet_times = (ui8*) calloc((size_t) 10000, sizeof(ui8));
    mef_channel_state_struct = (CHANNEL_STATE*) calloc((size_t) 1, sizeof(CHANNEL_STATE));
    
    // create MEF3 channel
    initialize_mef_channel_data(mef_channel_state_struct,
                                seconds_per_block,           // seconds per block
                                "sine-test"  , // channel name
                                0,// bit shift flag, set to 1 for neuralynx, to chop off 2 least-significant sample bits
                                0.0,           // low filt freq
                                9000.0,        // high filt freq
                                -1.0,           // notch filt freq
                                60.0,          // AC line freq
                                1.0,           // units conversion factor
                                "not entered ",// chan description
                                sampling_frequency, // starter freq for channel, make it as high or higher than actual freq to allocate buffers
                                1 * 1000000, // block interval, needs to be correct, this value is used for all channels
                                1,             // chan number
                                dir_name,      // absolute path of session
                                -6.0,                  // GMT offset
                                "not entered2",        // session description
                                "anon",                // anonymized subject name
                                "Mickey",                // subject first name
                                "Mouse",                 // subject second name
                                "",               // subject ID
                                "",           // institution
                                NULL,                  // level 1 password (technical data)
                                NULL,                  // level 2 password (subject data), must also specify level 1 password if specifying level 2
                                                       // level 1 and level 2 passwords should be different, if both are specified.
                                "not entered",        // study comments
                                "not entered",         // channel comments
                                0                      // secs per segment, 0 means no limit to segment size
                                );
    

    sine_amplitude = 20000.0;
    sine_frequency = 10.0;
    base_timestamp = 946684800000000;  // midnight, 1 January 2000
    
    // generate 10000 samples of sine wave
    for (i=0; i < 10000; i++)
    {
        samps[i] = (si4)(sine_amplitude * sin(2 * M_PI * i * sine_frequency / sampling_frequency));
        packet_times[i] = base_timestamp + (i * ((1e6)/sampling_frequency));  // extrapolate, putting into microseconds
    }
    
    // add buffered data to MEF channel
    // note: write_mef_channel_data can be called many times sequentially, so data can be compressed to MEF as the data
    // arrives.  The only caveat is that it must be pre-sorted in increasing time order.  
    // write_mef_channel_data() will not do any time sorting.
    write_mef_channel_data(mef_channel_state_struct, packet_times, samps, 10000, seconds_per_block, sampling_frequency);
    
    // all done, close MEF channel
    close_mef_channel(mef_channel_state_struct);
    
    /********************************  RECORDS ***************************************/
    
    // The following code is a simple demonstration of writing Note records.
    //
    // After the close_annotation() command, create_or_append_annotations() could be called again
    // and new records could then be appended to the same records files.
    
    // allocate struct
    annotation_state_struct = (ANNOTATION_STATE*) calloc((size_t) 1, sizeof(ANNOTATION_STATE));
    
    // create records (annotations) files
    create_or_append_annotations(annotation_state_struct, dir_name, -6.0, "not_entered");

    // manually write two "Note" type records
    write_annotation(annotation_state_struct, 946684800000000, "Note", "This is the text of the first note.");
    write_annotation(annotation_state_struct, 946684801000000, "Note", "This is the text of the second note.");
    
    // create a cursor record and write it to file.
    cursor_pointer = (MEFREC_Curs_1_0*) calloc((size_t) 1, sizeof(MEFREC_Curs_1_0));
    cursor_pointer->id_number = 1;
    cursor_pointer->trace_timestamp = 946684800000000;
    cursor_pointer->latency = 2000000;
    cursor_pointer->value = 10;
    sprintf(cursor_pointer->name, "My cursor");
    write_annotation(annotation_state_struct, 946684802000000, "Curs", cursor_pointer);
    
    // create an epoch record and write it to file.
    epoch_pointer = (MEFREC_Epoc_1_0*) calloc((size_t) 1, sizeof(MEFREC_Epoc_1_0));
    epoch_pointer->id_number = 1;
    epoch_pointer->timestamp = 946684803000000;
    epoch_pointer->end_timestamp = 946684804000000;
    epoch_pointer->duration = epoch_pointer->end_timestamp - epoch_pointer->timestamp;
    sprintf(epoch_pointer->epoch_type, "Generic");
    sprintf(epoch_pointer->text, "My example epoch");
    write_annotation(annotation_state_struct, 946684803000000, "Epoc", epoch_pointer);

    // close records files
    close_annotation(annotation_state_struct);
    
    // free allocated data
    free (annotation_state_struct);
    free (cursor_pointer);
    free (epoch_pointer);
    
    // Test reading the annotations we just wrote, by reading them and displaying them
    sprintf(temp_name, "%s.mefd/sine_test.rdat", dir_name);
    records_fps = read_MEF_file(NULL, temp_name, NULL, NULL, NULL, USE_GLOBAL_BEHAVIOR);
    if (records_fps->fp == NULL) {
        records_fps->fp = fopen(records_fps->full_file_name, "rb");
#ifndef _WIN32
        records_fps->fd = fileno(records_fps->fp);
#else
        records_fps->fd = _fileno(records_fps->fp);
#endif
    }
    show_records(records_fps);
    
    /********************************  END OF RECORDS ********************************/
    
    // free allocated data
    free (samps);
    free (packet_times);
    free (mef_channel_state_struct);
    
    return 1;
}
