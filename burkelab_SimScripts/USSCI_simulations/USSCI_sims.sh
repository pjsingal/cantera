# ./burkelab_SimScripts/USSCI_simulations/USSCI_sims.sh
fsz=8
fszxtick=7
fszytick=7
fszaxlab=9
lw=0.7
mw=0.5
msz=3.2
lgdw=1.0
lgdfsz=7

# python burkelab_SimScripts/USSCI_simulations/simulateshocktubeShao_generalLMRR.py \
# --figwidth 2.5 --figheight 2.5 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 3 &
# python burkelab_SimScripts/USSCI_simulations/simulateIDT_Shao_generalLMRR.py \
# --figwidth 2.5 --figheight 6.66667 --fsz $fszaxlab --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 7 --gridsz 9 --LMRtest 1 &
python burkelab_SimScripts/USSCI_simulations/simulateJSR_NH3_generalLMRR.py \
--figwidth 10 --figheight 5 --fsz $fszxtick --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
--lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 5.8 --gridsz 50 --date 'Oct22'
# python burkelab_SimScripts/USSCI_simulations/simulateJSR_H2O_generalLMRR.py \
# --figwidth 10 --figheight 5 --fsz $fszxtick --fszxtick $fszxtick --fszytick $fszytick --fszaxlab $fszaxlab \
# --lw $lw --mw $mw --msz $msz --lgdw $lgdw --lgdfsz 5.8 --gridsz 50 --date 'Oct22'
