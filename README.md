# KlaudSpeakerSimulator

The purpose of this repo is to try and recreate software like VituixCAD on Linux for simulating speakers.

All you need are the ZMA, FRD, and an idea of how your crossover circuits will look.

See TODO.md to see where I am on development.

All ZMA and FRD files are from: https://www.rjbaudio.com/Audiofiles/Driver%20FRD%20files.html
zma files: impedence data
	- columns: Frequency (Hertz), Impedence (Ohms), Phase (degrees)
frd files: sound pressure (loudness) level (SPL) data
	- columns: Frequency (Hertz), SPL (dB), Phase (degrees)
