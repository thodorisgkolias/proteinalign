DataExtract = function(data, atom = 'CA') {
     ## data must be an output from readpdb.
     ## selects the amino acid sequence and the coordinates columns.
     temp = data[data[, 3] == atom, c(4, 7, 8, 9)]
     return(temp)
}


SplitStr1 = function(st) {
     temp = unlist(strsplit(st, ' '))
     temp2 = temp[temp != '']
     return(temp2)	
}

AminoTransf = function(data) {
     ## Transforms the names of the amino acid sequence.
     
     freqs1 = c('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
                'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
                'TYR', 'VAL')
     freqs2 = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
                'F', 'P', 'S', 'T', 'W', 'Y', 'V')
     freqs = cbind(freqs1, freqs2)
     for (i in 1:length(levels(data$amino))) {
          aa = levels(data$amino)[i]
          rows = which(levels(data$amino)[i] == freqs[, 1])
          levels(data$amino) = sub(paste('^', aa, '$', sep =''),
                                   as.character(freqs[rows, 2]), 
                                   levels(data$amino))
     }
     return(data)
}

LoadPDB = function(pdb, atom = 'CA', chain = 'A') {
     if (nchar(pdb) != 4) {
          stop('No PDB file found')
     }
     ind1 = c(1, 7, 13, 17, 18, 22, 23, 27, 31, 39, 47, 55, 61, 73, 77, 79)
     ind2 = c(6, 11, 16, 17, 21, 22, 26, 27, 38, 46, 54, 60, 66, 76, 78, 80)
     file = paste('http://www.rcsb.org/pdb/files/', pdb, '.pdb', sep = '')
     data1 = readLines(file, n =-1)
     choose = substring(data1, 1, 6)
     atoms = data1[choose == 'ATOM  ']
     clean_data = lapply(atoms, substring, ind1, ind2)
     newd = lapply(clean_data, SplitStr1)
     rm_col = which(lapply(newd, length) != 12)
     if(length(rm_col) != 0) {
          for (i in rm_col) {
               newd[[i]] = newd[[i]][-7]
          }
     }
     newmat = matrix(unlist(newd), ncol = 12, byrow = TRUE)
     newmat = newmat[newmat[, 3] == atom, ]
     if (newmat[1, 5] != 'A') {
          chain = newmat[1, 5]
     }
     newmat = newmat[newmat[, 5] == chain, ]
     colnames(newmat) = c('type', 'seq', 'residue', 'amino', 'chain', 'res_no',
                          'x', 'y', 'z', 'occup', 'temp', 'elem_name')
     out = data.frame(newmat)
     out[, 7] = as.numeric(levels(out[, 7]))[out[, 7]]
     out[, 8] = as.numeric(levels(out[, 8]))[out[, 8]]
     out[, 9] = as.numeric(levels(out[, 9]))[out[, 9]]
     cat('', SplitStr1(data1[choose == 'HEADER'])[-1], '\n')
     
     out = DataExtract(data = AminoTransf(data = out),  atom = atom)
     return(out)
}



ReadPDB = function(dir){ 
     ## Read pdb format data from a directory
     
     ind1 = c(1, 7, 13, 17, 18, 22, 23, 27, 31, 39, 47, 55, 61, 73, 77, 79)
     ind2 = c(6, 11, 16, 17, 21, 22, 26, 27, 38, 46, 54, 60, 66, 76, 78, 80)
     data1 = readLines(dir, n =-1)
     choose = substring(data1, 1, 6)
     atoms = data1[choose == 'ATOM  ']
     clean_data = lapply(atoms, substring, ind1, ind2)
     newd = lapply(clean_data, SplitStr1)
     ll = lapply(newd, length)
     ll = unlist(ll)
     for (jj in 1:length(ll)){
          newd[[jj]] = newd[[jj]][1:min(ll)]
     }
     newmat = matrix(unlist(newd), ncol = ll, byrow = TRUE)
     out = data.frame(newmat)
     out[, 6] = as.numeric(levels(out[, 6]))[out[, 6]]
     out[, 7] = as.numeric(levels(out[, 7]))[out[, 7]]
     out[, 8] = as.numeric(levels(out[, 8]))[out[, 8]]
     out = AminoTransf(out[out[, 3] == 'CA',c(4, 6, 7, 8)])
     return(out)
}


PamTransf = function(pam) {
     res = round(10^(pam/10), 6)
     return(res)
}
