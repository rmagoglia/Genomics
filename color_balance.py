"""This script takes a list of adapter sequences for multiplexing
and prints out the color balance (green/red) ratio for each position
of the adapter. Ideal color balanced pools will have a green/red
ratio of 1.0 at each position. If any position has a value
of zero or inifinity, this pooling scheme cannot be used!

Usage: python3 color_balance.py"""



# Illumina adapters
D701 = "ATTACTCG"
D702 = "TCCGGAGA"
D703 = "CGCTCATT"
D704 = "GAGATTCC"
D705 = "ATTCAGAA"
D706 = "GAATTCGT"
D707 = "CTGAAGCT"
D708 = "TAATGCGC"
D709 = "CGGCTATG"
D710 = "TCCGCGAA"
D711 = "TCTCGCGC"
D712 = "AGCGATAG"

D501 = "TATAGCCT"
D502 = "ATAGAGGC"
D503 = "CCTATCCT"
D504 = "GGCTCTGA"
D505 = "AGGCGAAG"
D506 = "TAATCTTA"
D507 = "CAGGACGT"
D508 = "GTACTGAC"

# Read in list of adapters
# adapters = [D707,D709,D703,D704,D705,D706,D708,D710,D711,D712]
adapters = input("Please enter your adapter sequences, separated by a space: ")
adapters = adapters.split()

def check_length(adapters):
	lengths = []
	for i in adapters:
		lengths.append(len(i))

	lengths_set = set(lengths)

	if len(lengths_set) != 1:
		print("Error: Adapters have differing lengths.")
		return(False)
	else:
		return(True)

def check_colors(adapters):
	for i in range(0, len(adapters[0])):
		green = 0
		red = 0
		for j in adapters:
			if j[i] == "C" or j[i] == "A":
				red += 1
			elif j[i] == "G" or j[i] == "T":
				green +=1
			else:
				print("Error: Adapter \"" + j + "\" includes a noncanonical base.")
				return(False)

		try:
			ratio = str(round(green/red,1))
		except ZeroDivisionError:
			ratio = "inf"

		if ratio != "inf" and ratio != "0.0":
			status = "PASS"
		else:
			status = "FAIL"
		
		print("Position " + str(i+1) + ": green/red ratio of " + ratio + "     Status: " + status)
	return(True)

if check_length(adapters):
	check_colors(adapters)
