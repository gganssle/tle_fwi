#Velocity modeling toolkit

The scripts in this folder generate velocity models based on
explicit input. To change the model, simply alter the parameters
in the 'build model' commented section of the code.

Using alternative velocity model builders is perfectly
permissable for the included FWI algorithm. But keep in mind
that the wavelet envelope comparison algorithm included for
the tutorial is not sophisticated enough to handle wavelet
tuning, so use sparse vel models.
