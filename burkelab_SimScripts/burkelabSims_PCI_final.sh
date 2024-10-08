fsz=8
fszxtick=7
fszytick=7
fszaxlab=9
lw=0.7
mw=0.5
msz=3.2
lgdw=1.0
lgdfsz=7

# PCI PLOTS
# python burkelab_SimScripts/simulateIDT_Shao_4x1_PCI.py \
# --figwidth 2.5 --figheight 6.66667 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 9 &

# python burkelab_SimScripts/simulateJSR_H2O_3x1_PCI.py \
# --figwidth 2.5 --figheight 4.6 --fsz $fszxtick --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 50 &

# python burkelab_SimScripts/simulateJSR_NH3_3x1_PCI.py \
# --figwidth 2.5 --figheight 4.6 --fsz $fsz --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 50 &

# python burkelab_SimScripts/simulateshocktubeShao_1x1_PCI.py \
# --figwidth 2.5 --figheight 2.25 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 3 &

python burkelab_SimScripts/simulateflamespeedBurke_FromData_PCI.py \
--figwidth 2.5 --figheight 2.25 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --date 'Aug23' --slopeVal 0.05 --curveVal 0.05 &

python burkelab_SimScripts/simulateflamespeedRonney_0.6NH3_0.4H2_FromData.py \
--figwidth 2.5 --figheight 4.65 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz $lgdfsz --date 'Sep12' --slopeVal 0.05 --curveVal 0.05 --paper 'PCI' &

## USEFUL COMMANDS
# # Convert a cti to a yaml
# python interfaces\\cython\\cantera\\cti2yaml.py "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\09 Nitrogen\\Shrestha\\shrestha2018.cti" 

# # To make this file executable:
# chmod +x burkelab_SimScripts/burkelabSims.sh

# # To run this file:
# ./burkelab_SimScripts/burkelabSims_PCI_final.sh
