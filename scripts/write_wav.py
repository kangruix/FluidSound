# (c) 2024 Kangrui Xue
#
# write_wav.py
# Script for reading simulation output and writing to .wav file

import sys, os
import numpy as np, matplotlib.pyplot as plt
import soundfile as sf


def write_wav(filepath, samplerate):
    filename, ext = os.path.splitext(filepath)
    if ext == ".txt":
        data = np.loadtxt(filepath)
    elif ext == ".bin":
        data = np.fromfile(filepath, dtype=np.float32)
    else:
        raise Exception("Unsupported file extention: " + ext)

    data /= np.max(np.abs(data))  # normalize audio
    print(len(data) / samplerate, "s", "| Max:", np.max(data), "| Min:", np.min(data), "| Mean:", np.mean(data))

    wav_out = filename + ".wav"
    sf.write(wav_out, data, samplerate)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        write_wav(sys.argv[1], int(sys.argv[2]))
    else:
        print("Usage: write_wav.py <path to sim. output> <samplerate>")