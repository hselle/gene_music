from wavebender import *
print("=+=+=+=+=+=+=+=+=+=+WELCOME+=+=+=+=+=+=+=+=+=+=")

sample_dna2 = "ATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAATAAACG"
sample_dna3 = "ATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAAATTGCTGCCGCAGCGGGGCCGACGTCGAGTCAACAGAAGAAA"
fundemental = 440.0
def aminofy(dna):
    print('Commencing Sonification.........')
    aminos = list()
    for i in range(len(dna)):
        if i + 2 < len(dna):
            aminos.append(dna[i:i+3])
            i=i+3
    print('DNA aminofied.............')
    return aminos
print(aminofy(sample_dna2))

def sonify_aminos(amino_list):
    print('Sonifying Aminos.........')
    amplitude = .01
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

    pitch_dict = {
        "I": (square_wave(fundemental, amplitude=0.76*amplitude, ),),
        'L': (square_wave(fundemental * 2, amplitude=0.76*amplitude, ),),
        'V': (square_wave(fundemental * 1.5, amplitude=0.76*amplitude, ),),
        'F': (square_wave(fundemental * .5 , amplitude=0.76*amplitude, ),),
        'M': (square_wave(fundemental * 1.33, amplitude=0.76*amplitude, ),),
        "C": (square_wave(fundemental * 1.75, amplitude=0.76*amplitude, ),),
        'A': (square_wave(fundemental * 1.25, amplitude=0.76*amplitude, ),),
        'G': (square_wave(fundemental * .25, amplitude=0.76*amplitude, ),),
        'P': (square_wave(fundemental * .66, amplitude=0.76*amplitude, ),),
        'T': (square_wave(fundemental * .113, amplitude=0.76*amplitude, ),),
        "S": (square_wave(fundemental * 3, amplitude=0.76*amplitude, ),),
        'Y': (square_wave(fundemental * 3.5, amplitude=0.76*amplitude, ),),
        'W': (square_wave(fundemental * .25, amplitude=0.76*amplitude, ),),
        'Q': (square_wave(fundemental * 2.666, amplitude=0.76*amplitude, ),),
        'N': (square_wave(fundemental * 4, amplitude=0.76*amplitude, ),),
        "H": (square_wave(fundemental * 3.5, amplitude=0.76*amplitude, ),),
        'E': (square_wave(fundemental * .1, amplitude=0.76*amplitude, ),),
        'D': (square_wave(fundemental * 3.33, amplitude=0.76*amplitude, ),),
        'K': (square_wave(fundemental * 6.66, amplitude=0.76*amplitude, ),),
        'R': (square_wave(fundemental * 1, amplitude=0.76*amplitude, ),)

    }
    print("Made Pitch Dictionary...")
    chord = ()


    for amino in amino_list:
        if amino_chart[amino] == "STOP_CODON":
            ##stop
            print("-----Found Stop Codon------")
            return chord
        chord = chord + pitch_dict[amino_chart[amino]]
    return chord



chord = sonify_aminos(aminofy(sample_dna2))
channels = (chord,)

samples = compute_samples(channels, 44100 * 60 * 1)
write_wavefile("test_file2.wav", samples, 44100 * 60 * 1, nchannels=1)
