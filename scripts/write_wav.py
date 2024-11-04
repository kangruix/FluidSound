# (c) 2024 Kangrui Xue
#
# write_wav.py
# Script for reading simulation output and writing to .wav file

import sys
import numpy as np, matplotlib.pyplot as plt
import soundfile as sf


def write_wav(filename, samplerate, delim=" "):
    ftype = filename[-4:]
    if ftype == ".txt":
        data = np.loadtxt(filename)
    elif ftype == ".bin":
        with open(filename, "rb") as file:
            data = np.fromfile(file, dtype=np.float32)

    data /= np.max(np.abs(data)) * 1.01  # normalize audio
    data = (data * 32768).astype(np.int16)
    print(len(data) / samplerate, "s", "| Max:", np.max(data), "| Min:", np.min(data), "| Mean:", np.mean(data))

    out_file = filename.split(".")[0] + ".wav"
    sf.write(out_file, data, samplerate)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        write_wav(sys.argv[1], int(sys.argv[2]))
    else:
        print("Usage: write_wav.py <output file> <samplerate>")