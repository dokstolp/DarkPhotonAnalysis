The phoJetAnalyzer applies kinematic restrictions to the data, backgrounds, and signal samples. 

The analyzer requires:
1. A list of files from the dataset, an input value will determine which files to run over.
2. Histograms of the pileup vertices for the input dataset.
3. Scale factor 2d plot.
4. Signal k-factor plot.

To run:
./phoJetAnalyzer SampleName isMonteCarlo(1 for yes, 0 for no) iterator 

The iterator determines which files to run over as determined by line 771 of phoJetAnalyzer.h. By running all values of the iterator, one can analyze the entire dataset. The anlayzer breaks the dataset up into 20 files.
