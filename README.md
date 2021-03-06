# mefwriter
Code to simplify compressing data into the MEF 3.0 format.

This module of code (write_mef_channel files) is designed to be used in conjunction with the base
MEF 3.0 [library API](https://github.com/msel-source/meflib) (meflib and mefrec files).

Each channel of data, typically an electrode pair or voltage measurement, is a time-series array of data.

This module is not intended to be the definitive way to generate MEF 3.0 files.  It is simply one approach
that is generalized and is suitable for many applications.  The intention is to make MEF 3.0 usable,
without requiring a deep understanding of the details of the format.

Some sample data created using this module can be found [here](https://github.com/msel-source/sampledata).

General usage instructions:

This module assumes the timestamps are in increasing order, so if the original data isn't ordered properly,
it is the user's job to sort the data before using this module.  The threshold for discontinuity is specified
in write_mef_channel.h.  One tenth of a second, or 0.1 seconds, is the default.

The example program, named sine-test, is provided to demonstrate ease of use of this module.  The basic
pattern is 1) a channel is created, 2) data is added to the channel, and 3) the channel is closed.  The data
need not be added all at once; step 2 in the above list can be done many times, in between steps 1 and 3.  Though 
not shown in the example program, more data can be added to a pre-existing (and closed) channel using 
the append function, and this will create a new segment of data in the channel.

The example program is in both C and C#.  With C# things are a little more tricky, since the base MEF 3.0
API and the write_mef_channel module need to be compiled in C, and exported as a .dll.  The MSEL lab
isn't officially supporiting C#, but the code is provided to show an example of use.

Compiling this code under C, which is extensively tested on Mac OS X, is relatively simple.  The files
in this repository, as well as the files in the meflib repository, can be compiled using "gcc *.c".

Detailed usage notes:

Multiple channels can be created within one program.  The example programs only show one channel, but
the code is fully object-oriented in the sense that many channels can be created in one program.  While we don't
guarantee thread safety (such as when adding data to multiple channels simultaneously), the only issue we are 
aware of is the timestamp-offset generation issue, where every MEF 3.0 recording session has one (and only one)
offset for timestamp encryption.  Generally we implement timestamp offsetting even in cases where encryption 
is not used.  The offset is generated in new channels in the process_filled_block() subroutine in 
write_mef_channel.c, which is where a mutex could be added.

Do not add data to the same channel simultaneously from multiple threads.  There is no good reason to do that
anyway, since data might not be ordered properly.

When designing RED (range encoded differences) block sizes for channels, our research suggests optimal 
compression occurs in the 20000 to 30000 samples-per-block range.  Of course the nature of the data
will determine how well the compression works.  We like to use 1 second blocks for high frequency data,
such as 30kHz.  For low frequency data, such as 250 or 500 Hz, we like 15 second blocks, which seems to
be a good tradeoff between compression and ease-of-use for decompressing data later on.  If blocks get
too large, then a lot of data has to be decompressed to find what you are looking for.

This software is licensed under the Apache software license 2.0. See [LICENSE](./LICENSE) for details.
