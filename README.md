# BasicGDIC

basicgdic.m is a matlab function that provides a graphical user interface to perform 2D Digital Image Correlation (DIC). The goal in DIC is to find the displacement field that occurred between two images of the same object that are taken at different times. Typically, the object contains some form of speckle or texture to help this procedure. This method attempts to minimize the image residuals globally over the entire region of interest, which is different from the local formulation normally found in most commercial DIC tools. The tool is called basic since it was never designed to be a full-fledged DIC software suit. The main reason for developing this tool was for educational purposes. Admittedly, the amount of features and options available in this tool counteract the name, but this is all subjective.

Being developed in matlab, with the focus on educational purposes has a few consequences:
- The tool is very open, not only open-source, but also the data and algorithms are designed to be easily accessible.
- The tool has lots of options that are commonly hidden away, and are probably superfluous for normal DIC users.
- The tool is very honest about the quality of the data, not hiding any residuals nor smoothing any results.
- The tool allows a DIC calculation to be initialized with a previous DIC result, this allows progressive refinements which is very useful for very difficult to converge cases, at the expense of time.
- The tool is NOT fast, and is waisting many CPU cycles on needlessly moving data in and out of memory.
- The tool requires lots of memory. Many data sources are kept in memory to simplify the data flow of the tool. Consequently, this tool requires a lot of available memory (e.g. 8 GB) even for moderate size images (e.g. 2 MPixels).

Finally, for more information, please start the GUI and find the HELP section.