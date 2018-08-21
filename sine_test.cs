// Multiscale Electrophysiology Format (MEF) version 3.0
// Copyright 2018, Mayo Foundation, Rochester MN. All rights reserved.
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

// This code is a simple example in C# of how to create a MEF 3.0 channel, add some data to it,
// then close the channel.  This example assumes the write_mef_channel c code has been compiled
// as a .dll called "DllTest1".

using System;
using System.Runtime.InteropServices;

namespace SineTestApp
{
    class Program
    {

        public unsafe struct CHANNEL_STATE
        {
            public fixed byte first[5000];   // CHANNEL_STATE size is 2576, according to sizeof() in c.
        }


#if WIN64
        [DllImport("DllTest1")]
#else
        [DllImport("DllTest1", CallingConvention = CallingConvention.Cdecl)]
#endif
        protected static extern int initialize_meflib_dll();

#if WIN64
        [DllImport("DllTest1")]
#else
        [DllImport("DllTest1", CallingConvention = CallingConvention.Cdecl)]
#endif
        protected static extern int close_mef_channel(ref CHANNEL_STATE inputs);

#if WIN64
        [DllImport("DllTest1")]
#else
        [DllImport("DllTest1", CallingConvention = CallingConvention.Cdecl)]
#endif
        protected static extern int update_metadata_dll(ref CHANNEL_STATE inputs);

#if WIN64
        [DllImport("DllTest1")]
#else
        [DllImport("DllTest1", CallingConvention = CallingConvention.Cdecl)]
#endif
        protected static extern int write_mef_channel_data(ref CHANNEL_STATE inputs,
                                                       ulong[] packet_times,
                                                       int[] samps,
                                                       ulong n_packets_to_process,
                                                       double secs_per_block,
                                                       double sampling_frequency);

#if WIN64
        [DllImport("DllTest1")]
#else
        [DllImport("DllTest1", CallingConvention = CallingConvention.Cdecl)]
#endif
        protected static extern int flush_mef_channel(ref CHANNEL_STATE inputs);

#if WIN64
        [DllImport("DllTest1")]
#else
        [DllImport("DllTest1", CallingConvention = CallingConvention.Cdecl)]
#endif
        protected static extern int initialize_mef_channel_data(ref CHANNEL_STATE inputs,
                                                            double secs_per_block,
                                                            string chan_map_name,
                                                            int bit_shift_flag,
                                                            double low_frequency_filter_setting,
                                                            double high_frequency_filter_setting,
                                                            double notch_filter_frequency,
                                                            double AC_line_frequency,
                                                            double units_conversion_factor,
                                                            string channel_description,
                                                            double sampling_frequency,
                                                            long block_interval,
                                                            int chan_num,
                                                            string mef3_session_directory,
                                                            float gmt_offset,
                                                            string session_description,
                                                            string anonymized_subject_name,
                                                            string subject_first_name,
                                                            string subject_second_name,
                                                            string subject_id,
                                                            string institution,
                                                            string mef_3_level_1_password,
                                                            string mef_3_level_2_password,
                                                            string study_comments,
                                                            string channel_comments,
                                                            ulong num_secs_per_segment
                                                            );

#if WIN64
        [DllImport("DllTest1")]
#else
        [DllImport("DllTest1", CallingConvention = CallingConvention.Cdecl)]
#endif
        protected static extern int append_mef_channel_data(ref CHANNEL_STATE inputs,
        string chan_map_name,
        int new_segment_number,
        string mef_3_level_1_password,
        string mef_3_level_2_password,
        string mef3_session_directory,
        ulong num_secs_per_segment,
        int bit_shift_flag);

        static void Main(string[] args)
        {

            int i;
            double sampling_frequency;
            double sine_amplitude;
            double sine_frequency;
            double seconds_per_block;
            CHANNEL_STATE mef_channel_state_struct;
            int[] samps;
            ulong[] packet_times;
            ulong base_timestamp;
            string dir_name;
            string chan_name;

            // initialize variables
            sampling_frequency = 1000.0;   // Hz
            seconds_per_block = 1.0;
            dir_name = "c:/sinetest2";
            chan_name = "sinetest2";

            mef_channel_state_struct = new CHANNEL_STATE();

            // allocate buffers
            samps = new int[10000];
            packet_times = new ulong[10000];

            // set base data
            sine_amplitude = 20000.0;
            sine_frequency = 10.0;
            base_timestamp = 946684800000000;  // midnight, 1 January 2000, in microseconds

            // initialize MEF3 library
            initialize_meflib_dll();

            Console.WriteLine("Done initializing MEF library.");

            // create MEF3 channel
            initialize_mef_channel_data(ref mef_channel_state_struct,
                                        seconds_per_block,           // seconds per block
                                        chan_name, // channel name
                                        0,// bit shift flag, set to 1 for neuralynx, to chop off 2 least-significant sample bits
                                        0.0,           // low filt freq
                                        9000.0,        // high filt freq
                                        -1.0,           // notch filt freq
                                        60.0,          // AC line freq
                                        1.0,           // units conversion factor
                                        "not entered", // chan description
                                        sampling_frequency, // starter freq for channel, make it as high or higher than actual freq to allocate buffers
                                        1 * 1000000, // block interval, needs to be correct, this value is used for all channels
                                        1,             // chan number
                                        dir_name,      // absolute path of session
                                        (float)-6.0,                  // GMT offset
                                        "not entered",        // session description
                                        "anon",                // anonymized subject name
                                        "First",                // subject first name
                                        "Last",                 // subject second name
                                        "",               // subject ID
                                        "",           // institution
                                        null,                  // level 1 password (technical data)
                                        null,                  // level 2 password (subject data), must also specify level 1 password if specifying level 2
                                        "not entered",        // study comments
                                        "not entered",         // channel comments
                                        0                      // secs per segment, 0 means no limit to segment size
                                        );

            Console.WriteLine("Done initializing MEF channels.");

            // generate 10000 samples of sine wave
            for (i = 0; i < 10000; i++)
            {
                samps[i] = (int)(sine_amplitude * Math.Sin(2 * Math.PI * i * sine_frequency / sampling_frequency));
                packet_times[i] = (ulong)(base_timestamp + (i * ((1e6) / sampling_frequency)));  // extrapolate, putting into microseconds
            }

            // add buffered data to MEF channel    
            write_mef_channel_data(ref mef_channel_state_struct, packet_times, samps, 10000, seconds_per_block, sampling_frequency);
  
            Console.WriteLine("Done adding data to MEF channel.");

            // close MEF channel
            close_mef_channel(ref mef_channel_state_struct);

            Console.WriteLine("Closed MEF channel.");
            Console.WriteLine("Press enter to end the program.");

            Console.ReadLine();
        }
    }
}
