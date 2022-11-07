APP_NAME      := analysis

SOURCE_FILES  := analysis.cc
USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := yes

include $(HADDIR)/hades.def.mk

HYDRA_LIBS    += -lDst


.PHONY:  default
default: clean build install

# any special lib needed?
# libDst is not in the default
HYDRA_LIBS    += -lMLP -lMinuit -lTreePlayer -lTMVA -lXMLIO -lm    #-lTMVAGui #needed for KZeroTrees
#HYDRA_LIBS    += -lEvtChara    #needed for embedding

include $(HADDIR)/hades.app.mk
