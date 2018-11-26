from wavebender import *
from itertools import chain, cycle
import time
import sys

#global progress = 0
def aminofy(dna):
    print('Commencing Sonification.........')
    aminos = list()
    for i in range(len(dna)):
        if i + 2 < len(dna):
            aminos.append(dna[i:i+3])
            i=i+3
    print('DNA aminofied.............')
    return aminos

def sonify_aminos(amino_list, fundemental, amplitude, FRAMERATE):
    print('Sonifying Aminos.........')

    amino_chart = {
        "ATT": "I",
        'ATC': 'I',
        'ATA': 'I',
        'CTT': 'L',
        'CTC': 'L',
        'CTA': 'L',
        'CTG': 'L',
        'TTA': 'L',
        'TTG': 'L',
        'GTT': 'V',
        'GTC': 'V',
        'GTA': 'V',
        'GTG': 'V',
        'TTT': 'F',
        'TTC': 'F',
        'ATG': 'M',
        'TGT': 'C',
        'TGC': 'C',
        'GCT': 'A',
        'GCC': 'A',
        'GCA': 'A',
        'GCG': 'A',
        'GGT': 'G',
        'GGC': 'G',
        'GGA': 'G',
        'GGG': 'G',
        'CCT': 'P',
        'CCC': 'P',
        'CCA': 'P',
        'CCG': 'P',
        'ACT': 'T',
        'ACC': 'T',
        'ACA': 'T',
        'ACG': 'T',
        'TCT': 'S',
        'TCC': 'S',
        'TCA': 'S',
        'TCG': 'S',
        'AGT': 'S',
        'AGC': 'S',
        'TAT': 'Y',
        'TAC': 'Y',
        'TGG': 'W',
        'CAA': 'Q',
        'CAG': 'Q',
        'AAT': 'N',
        'AAC': 'N',
        'CAT': 'H',
        'CAC': 'H',
        'GAA': 'E',
        'GAG': 'E',
        'GAT': 'D',
        'GAC': 'D',
        'AAA': 'K',
        'AAG': 'K',
        'CGT': 'R',
        'CGC': 'R',
        'CGA': 'R',
        'CGG': 'R',
        'AGA': 'R',
        'AGG': 'R',
        'TAA': 'STOP_CODON',
        'TAG': 'STOP_CODON',
        'TGA': 'STOP_CODON'
    }
    print("Creating waves......")
    iwave = (sine_wave(fundemental, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    lwave = (sine_wave(fundemental * 2, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    vwave = (sine_wave(fundemental * 1.5, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    fwave = (sine_wave(fundemental * .5 , amplitude=0.76*amplitude, framerate=FRAMERATE),)
    mwave = (sine_wave(fundemental * 1.33, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    cwave = (sine_wave(fundemental * 1.75, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    awave = (sine_wave(fundemental * 1.25, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    gwave = (sine_wave(fundemental * .25, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    pwave = (sine_wave(fundemental * .66, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    twave = (sine_wave(fundemental * .113, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    swave = (sine_wave(fundemental * 3, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    ywave = (sine_wave(fundemental * 3.5, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    wwave = (sine_wave(fundemental * .25, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    qwave = (sine_wave(fundemental * 2.666, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    nwave = (sine_wave(fundemental * 4, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    hwave = (sine_wave(fundemental * 3.5, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    ewave = (sine_wave(fundemental * .1, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    dwave = (sine_wave(fundemental * 3.33, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    kwave = (sine_wave(fundemental * 6.66, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    rwave = (sine_wave(fundemental * 1, amplitude=0.76*amplitude, framerate=FRAMERATE),)
    print("Waves Created")
    pitch_dict = {
        "I": iwave,
        'L': lwave,
        'V': vwave,
        'F': fwave,
        'M': mwave,
        "C": cwave,
        'A': awave,
        'G': gwave,
        'P': pwave,
        'T': twave,
        "S": swave,
        'Y': ywave,
        'W': wwave,
        'Q': qwave,
        'N': nwave,
        "H": hwave,
        'E': ewave,
        'D': dwave,
        'K': kwave,
        'R': rwave

    }
    print("Made Pitch Dictionary...")
    chord_list = list()
    chord = ()
    stop_codon_count = 0
    for amino in amino_list:

        if amino_chart[amino] == "STOP_CODON":
            stop_codon_count = stop_codon_count + 1
            print("-----Found Stop Codon------" + str(stop_codon_count))
            chord_list.append(chord)
            chord = ()
        else:
            chord = chord + pitch_dict[amino_chart[amino]]
    ##chord_list.append(chord)
    ##print("Here is the chord list: " + str(chord_list))
    return chain(chain((_chord) for _chord in chord_list))
def getProgress():
    print("=+=+=+=+=+=+=+=+=+=+WELCOME+=+=+=+=+=+=+=+=+=+=")
    global progress
    while progress < len(gene_string):
        print("WORKING-_-_-__-__-_-_-__-" + str(progress))
        time.sleep(.3)
def main():
    parser = argparse.ArgumentParser(prog="gene_audio")
    parser.add_argument('in_file', help="The file containing the genetic seqeunce")
    parser.add_argument('out_file', help="The generated wav file")
    args = parser.parse_args()

    timestart = time.clock()
    print("=+=+=+=+=+=+=+=+=+=+WELCOME+=+=+=+=+=+=+=+=+=+=")


    f = open(args.in_file, 'r')
    gene_string = f.read()[:-1]
    print("length of collo: " + str(len(gene_string)))
    gene_string = gene_string

    fundemental = 440.0
    amplitude = .001
    FRAMERATE = 22050

    chords = sonify_aminos(aminofy(gene_string), fundemental, amplitude, int(FRAMERATE * .1 * 1))
    channels = chords
    print("Here are the channels: " + str(channels))
    print("length of collo: " + str(len(gene_string)))
    samples = compute_samples(channels, int(FRAMERATE * .1 * 1))

    write_wavefile(args.out_file, samples, int(FRAMERATE * .1), framerate = FRAMERATE, nchannels=1, bufsize = 128)

    print(time.clock() - timestart)
if __name__ == '__main__':
    main()
