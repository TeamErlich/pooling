#!/usr/bin/sage

###################################
#
# This file implements the recursive argument described in the supplementary material: 
# "computing the average disjunctness of Reed-Solomon Codes."
#
##################################
q=16
n=6

def getq2(m,r,x):
	return binomial(x, 2)/m

def getq1(m,r,x):
	return (q*x + x*(q-x)) /m

def qs(m,r,x):
	q1 = getq1(m,r,x)
	q2 = getq2(m,r,x)
	return [1 - q1-q2, q1,q2]

def compute( m,r,x, table_so_far ):
	if r < x/2:
		ret = 0
	elif m==0:
		ret = 0 
	elif x==1 and r==1:
		ret = (q + (q-1))/m
	elif x==2 and r==1:
		ret = 1/m
	elif x==0:
		ret = 1
	else:
		q0,q1,q2 = qs(m,r,x)
		if (m-1,r-1,x) not in table_so_far.keys():
			compute(m-1,r-1,x, table_so_far)
		if (m-1,r-1,x-1) not in table_so_far.keys():
			compute(m-1,r-1,x-1, table_so_far)
		if (m-1,r-1,x-2) not in table_so_far.keys():
			compute(m-1,r-1,x-2, table_so_far)
		ret = q0 * table_so_far[(m-1,r-1,x)] + q1 * table_so_far[(m-1,r-1,x-1)] + q2* table_so_far[(m-1,r-1,x-2)]
	table_so_far[(m,r,x)] = ret
	return ret

def count(d,N):
	table = {}
	ret = compute( N, d, n, table)
	return ret

F = open("outfiles/theoretical.dat", "w")
for d in range(1, 15):
	F.writelines(str(d) + "\t")
	for N in [4*96,5*96,6*96]:
		val = 1 - N*float(count(d,N))
		F.writelines(str(val) + "\t")
	F.writelines("\n")
F.close()
