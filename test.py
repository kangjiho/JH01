class Namedef Name@indenter("Scheme Name")
def Scheme Name_indent(c, doc):
	sci = scintilla.Scintilla(doc)
	editor = scintilla.Scintilla(pn.CurrentDoc())():
	:
	
	raw_input()
	input()
	
	
    #{ POST TL
    if POST==True:

      # FILENAME
      Filename = FO_DIR+'/(%2d,%2d)_TL.csv'%(K,L)

      #{ HEADER
      try: fo = open(Filename,'r')
      except: # File doesn't exist -> This is the first time. Make header.
        temp  = 't, '
        temp += 'Epsil, F_TLk, dEpsil, F_TLc, F_TL, M_TL\n'
        fo = open(Filename,'w')
        fo.write(temp)
        fo.close()
      #}

      #{ F_TLk, F_TLc
      F_TLk = F_TLc = 0.
      if Epsilon > 0:
        F_TLk = +Kh *  Epsilon
        F_TLc = +Ch * DEpsilon
      #}

      #{ WRITE
      temp  = '%e, '%t
#              1   2   3   4   5   6    : 1        2      3         4      5      6
      temp += '%e, %e, %e, %e, %e, %e\n'%(Epsilon, F_TLk, DEpsilon, F_TLc, Force, Moment )
      fo = open(Filename,'a')
      fo.write(temp)
      fo.close()
      #}
    #}





