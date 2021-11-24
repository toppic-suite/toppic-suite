#!/bin/bash
start=`date +%s`

./test_graph -i mass_graph_mods.txt -k -p 0 -j 20 H4.fasta 071210_070610His0Gy070210H4_H061010A_msdeconv.msalign

end=`date +%s`
echo "Runnig time: " $((end-start)) " seconds"
