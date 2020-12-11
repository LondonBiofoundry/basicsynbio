import basicsynbio as bsb

mypart = bsb.BASIC_SEVA_PARTS['v0.1']["27"]

#print('\n')
#print(mypart)
#print('\n')

mylinker = bsb.BASIC_BIOLEGIO_LINKERS['v0.1']["LMP"]

#print('\n')
#print(mylinker)
#print('\n')

coll = bsb.BASIC_SEVA_PARTS['v0.1']

#print(coll)

myassem = bsb.BasicAssembly(
            "test",
            bsb.BASIC_SEVA_PARTS['v0.1']["27"],
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["UTR1-RBS1"],
            bsb.BASIC_SEVA_PARTS['v0.1']["19"],
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["UTR1-RBS2"],
        )
myassem.return_clip_reactions()
