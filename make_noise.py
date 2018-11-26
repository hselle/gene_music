from wavebender import *
import time
timestart = time.clock()
print("=+=+=+=+=+=+=+=+=+=+WELCOME+=+=+=+=+=+=+=+=+=+=")

sample_dna2 = "ATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAATAAACG"

file_name = 'collogen.txt'
f = open(file_name, r)
fundemental = 440.0
x = 5.50/len(collogen)
amplitude = .00005
FRAMERATE = 22050
def aminofy(dna):
    print('Commencing Sonification.........')
    aminos = list()
    for i in range(len(dna)):
        if i + 2 < len(dna):
            aminos.append(dna[i:i+3])
            i=i+3
    print('DNA aminofied.............')
    return aminos

def sonify_aminos(amino_list):
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
    iwave = (damped_wave(fundemental, amplitude=0.76*amplitude, length=FRAMERATE),)
    lwave = (damped_wave(fundemental * 2, amplitude=0.76*amplitude, length=FRAMERATE),)
    vwave = (damped_wave(fundemental * 1.5, amplitude=0.76*amplitude, length=FRAMERATE),)
    fwave = (damped_wave(fundemental * .5 , amplitude=0.76*amplitude, length=FRAMERATE),)
    mwave = (damped_wave(fundemental * 1.33, amplitude=0.76*amplitude, length=FRAMERATE),)
    cwave = (damped_wave(fundemental * 1.75, amplitude=0.76*amplitude, length=FRAMERATE),)
    awave = (damped_wave(fundemental * 1.25, amplitude=0.76*amplitude, length=FRAMERATE),)
    gwave = (damped_wave(fundemental * .25, amplitude=0.76*amplitude, length=FRAMERATE),)
    pwave = (damped_wave(fundemental * .66, amplitude=0.76*amplitude, length=FRAMERATE),)
    twave = (damped_wave(fundemental * .113, amplitude=0.76*amplitude, length=FRAMERATE),)
    swave = (damped_wave(fundemental * 3, amplitude=0.76*amplitude, length=FRAMERATE),)
    ywave = (damped_wave(fundemental * 3.5, amplitude=0.76*amplitude, length=FRAMERATE),)
    wwave = (damped_wave(fundemental * .25, amplitude=0.76*amplitude, length=FRAMERATE),)
    qwave = (damped_wave(fundemental * 2.666, amplitude=0.76*amplitude, length=FRAMERATE),)
    nwave = (damped_wave(fundemental * 4, amplitude=0.76*amplitude, length=FRAMERATE),)
    hwave = (damped_wave(fundemental * 3.5, amplitude=0.76*amplitude, length=FRAMERATE),)
    ewave = (damped_wave(fundemental * .1, amplitude=0.76*amplitude, length=FRAMERATE),)
    dwave = (damped_wave(fundemental * 3.33, amplitude=0.76*amplitude, length=FRAMERATE),)
    kwave = (damped_wave(fundemental * 6.66, amplitude=0.76*amplitude, length=FRAMERATE),)
    rwave = (damped_wave(fundemental * 1, amplitude=0.76*amplitude, length=FRAMERATE),)
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
    chord = ()

    stop_codon_count = 0
    for amino in amino_list:

        if amino_chart[amino] == "STOP_CODON":
            stop_codon_count = stop_codon_count + 1
            print("-----Found Stop Codon------" + str(stop_codon_count))
        else:
            chord = chord + pitch_dict[amino_chart[amino]]
    return chord

chord = sonify_aminos(aminofy(f.read()))
channels = (chord,)

samples = compute_samples(channels, FRAMERATE * 1 * 1)
write_wavefile("collogen.wav", samples, FRAMERATE * 1 * 1, framerate = FRAMERATE, nchannels=1, bufsize = 2048)
print(time.clock() - timestart)
