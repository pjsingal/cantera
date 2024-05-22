

python burkelab_SimScripts/simulateIDT_Shao_1x4_ESSCI.py
python burkelab_SimScripts/simulateIDT_Shao_4x1_PCI.py
python burkelab_SimScripts/simulateshocktubeShao_1x1_ESSCI.py
python burkelab_SimScripts/simulateshocktubeShao_1x1_PCI.py
python burkelab_SimScripts/simulateJSR_H2O_1x3_ESSCI.py
python burkelab_SimScripts/simulateJSR_H2O_3x1_PCI.py
python burkelab_SimScripts/simulateJSR_NH3_1x3_ESSCI.py
python burkelab_SimScripts/simulateJSR_NH3_3x1_PCI.py

"""
# To make this file executable:
chmod +x burkelab_SimScripts/burkelabSims.sh

# To run this file:
./burkelab_SimScripts/burkelabSims.sh

#1x4
# figsize = (165/25.4, 50/25.4) # (width, height)
# fig, ax = plt.subplots(1, 4, figsize=figsize)
figsize = (6.5, 1.5) # (width, height)
fig, ax = plt.subplots(1, 4, figsize=figsize)

#4x1
# figsize = (90/25.4, 110/25.4) # (width, height)
# fig, ax = plt.subplots(4, 1, figsize=figsize)
figsize = (3.5, 7.5) # (width, height)
fig, ax = plt.subplots(4, 1, figsize=figsize)

#1x3
# figsize = (165*0.75/25.4, 50/25.4) # (width, height)
# fig, ax = plt.subplots(1, 3, figsize=figsize)
# f, ax = plt.subplots(1, 3, figsize=(6, 1.5)) 
f, ax = plt.subplots(1, 3, figsize=(6, 1.5))

#3x1
# figsize = (90/25.4, 110*0.75/25.4) # (width, height)
# fig, ax = plt.subplots(3, 1, figsize=figsize)
# f, ax = plt.subplots(3, 1, figsize=(9, 5))
f, ax = plt.subplots(3, 1, figsize=(3.5, 7))

#1x1
# figsize = (90/25.4, 50/25.4) # (width, height)
# fig, ax = plt.subplots(1, 1, figsize=figsize)
# fig, ax = plt.subplots(1, 1, figsize=(9, 5))
fig, ax = plt.subplots(1, 1, figsize=(3.8, 2))


lw=0.7     #linewidth
mw=0.5     #marker edge width
msz=3.5    #marker size
dpi=1000
lgdw=0.6   #legend width
lgdfsz=7   #legend font size

"""