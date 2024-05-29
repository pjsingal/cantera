fsz=8
fszxtick=7
fszytick=7
fszaxlab=8
lw=0.7
mw=0.5
msz=2.5
lgdw=1.0
lgdfsz=7

python burkelab_SimScripts/simulateIDT_Shao_1x4_ESSCI.py \
--figwidth 6 --figheight 1.66667 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 9 &

python burkelab_SimScripts/simulateIDT_Shao_4x1_PCI.py \
--figwidth 3.5 --figheight 6.66667 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 9 &

python burkelab_SimScripts/simulateJSR_H2O_1x3_ESSCI.py \
--figwidth 6 --figheight 1.66667 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 5 --gridsz 50 &

python burkelab_SimScripts/simulateJSR_H2O_3x1_PCI.py \
--figwidth 3.5 --figheight 5 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 50 &

python burkelab_SimScripts/simulateJSR_NH3_1x3_ESSCI.py \
--figwidth 6 --figheight 1.66667 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 50 &

python burkelab_SimScripts/simulateJSR_NH3_3x1_PCI.py \
--figwidth 3.5 --figheight 5 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 50 &

python burkelab_SimScripts/simulateshocktubeShao_1x1_ESSCI.py \
--figwidth 3.5 --figheight 1.66667 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz 2.5 --lgdw $lgdw --lgdfsz $lgdfsz --gridsz 14 &

python burkelab_SimScripts/simulateshocktubeShao_1x1_PCI.py \
--figwidth 3.5 --figheight 1.66667 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz 2.5 --lgdw $lgdw --lgdfsz $lgdfsz --gridsz 14




# # To make this file executable:
# chmod +x burkelab_SimScripts/burkelabSims.sh

# # To run this file:
# ./burkelab_SimScripts/burkelabSims.sh