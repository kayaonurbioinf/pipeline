 Have the pipeline.py in the same folder as you unpacked the MHC I and MHC II folders to. Create two additional folders called fas_out and my_annotations. Your parsed annotations (.json) file will be created in my_annotations, domains files will be generated in my_featuretypes and the MHC result tables will be generated in the same folder as the script. You must also have the input file in the same folder as the pipeline. In total you must have the following in the same directory:

    my_annotations (folder)
    my_featuretypes (folder)
    MHC I and MHC II folders (you can just unpack the tar.gz files into the same folder as above)
    your input file (fasta)

You can use the command “predict” to do MHC binding predictions and run FAS on the file. The command must be in the following form if you only need MHC I:

    python pipeline.py predict my_input.fasta mhc_i 10

For only MHC II:

    python pipeline.py predict my_input.fasta mhc_ii 15

For both:

    python pipeline.py predict my_input.fasta both 10 15

When you enter “both” as the tool argument, you will need to enter 2 lengths afterwards. The first length will be passed to MHC I and the second to MHC II. If you want to look for a range of lengths you can enter the length this way:

    python pipeline.py predict my_input.fasta mhc_i 10  12-14

When you enter a length in the x-y form, the pipeline will look for all lengths including x and y. E.g. if you enter 9-14, the pipeline will look for binding sites of lenghts 9, 10, 11, 12, 13 and 14.

Additionally, you can use the command “info” for a number of auxiliary purposes:

    python pipeline.py info predict

will display a text on how exactly the command predict is supposed to be used:

python pipeline.py ls will print a list of files in the current working directory. Note that the pipeline extends its working directory by a few directories to import certain modules. This might be useful of you are having trouble with getting the pipeline import a module or if for some reason the pipeline won't see your input file.

    python pipeline.py cwd

will display the current working directory.
