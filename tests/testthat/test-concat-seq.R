test_that("Test Case 1: 9-mer HCS (single host), with CSV output",{
    HCS_1host<- concat_conserved_kmer(proteins_1host)

    expect_equal(nrow(HCS_1host),27)
    expect_equal(ncol(HCS_1host),3)

    col_names<-colnames(HCS_1host)
    expect_equal(col_names[1],'HCS')
    expect_equal(col_names[2],'Position')
    expect_equal(col_names[3],'Sequence')

    #randomly pick 3 rows to check
    #row 1
    firstRow <- HCS_1host[1,]
    expect_equal(firstRow$HCS, 'HCS_A_1')
    expect_equal(firstRow$Position, '1-18')
    expect_equal(firstRow$Sequence, 'MSTNPKPQRKTKRNTNRR')

    #row 10
    tenthRow <- HCS_1host[10,]
    expect_equal(tenthRow$HCS, 'HCS_B_2')
    expect_equal(tenthRow$Position, '154-166')
    expect_equal(tenthRow$Sequence, 'FRAAVCTRGVAKA')

    #row 27
    row27th <- HCS_1host[27,]
    expect_equal(row27th$HCS, 'HCS_B_19')
    expect_equal(row27th$Position, '587-603')
    expect_equal(row27th$Sequence, 'RLKPTLHGPTPLLYRLG')
})

test_that("Test Case 2: 9-mer CCS (single host), with CSV output",{
    CCS_1host<- concat_conserved_kmer(proteins_1host, conservationLevel = 'CCS')

    expect_equal(nrow(CCS_1host),1)
    expect_equal(ncol(CCS_1host),3)

    col_names<-colnames(CCS_1host)
    expect_equal(col_names[1],'CCS')
    expect_equal(col_names[2],'Position')
    expect_equal(col_names[3],'Sequence')

    #row 1
    firstRow <- CCS_1host[1,]
    expect_equal(firstRow$CCS, 'CCS_A_1')
    expect_equal(firstRow$Position, '1-18')
    expect_equal(firstRow$Sequence, 'MSTNPKPQRKTKRNTNRR')

})

test_that("Test Case 3: 9-mer HCS (single host), with FASTA output",{
    HCS_1host<- concat_conserved_kmer(proteins_1host, output_type = "fasta")

    expect_equal(nrow(HCS_1host),54)
    expect_equal(ncol(HCS_1host),1)

    #randomly pick 3 FASTA to check
    #FASTA 1
    expect_equal(HCS_1host[1,1], '>HCS_A_1')
    expect_equal(HCS_1host[2,1], 'MSTNPKPQRKTKRNTNRR')

    #FASTA 10
    expect_equal(HCS_1host[19,1], '>HCS_B_2')
    expect_equal(HCS_1host[20,1], 'FRAAVCTRGVAKA')

    #FASTA 27
    expect_equal(HCS_1host[53,1], '>HCS_B_19')
    expect_equal(HCS_1host[54,1], 'RLKPTLHGPTPLLYRLG')
})

test_that("Test Case 4: 9-mer CCS (single host), with FASTA output",{
    CCS_1host<- concat_conserved_kmer(proteins_1host, conservationLevel = 'CCS', output_type = "fasta")

    expect_equal(nrow(CCS_1host),2)
    expect_equal(ncol(CCS_1host),1)

    expect_equal(CCS_1host[1,1], '>CCS_A_1')
    expect_equal(CCS_1host[2,1], 'MSTNPKPQRKTKRNTNRR')

})