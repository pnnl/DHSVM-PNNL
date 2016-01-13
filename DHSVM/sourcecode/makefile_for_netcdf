# 
# SUMMARY:      makefile for DHSVM
# USAGE:        make DHSVM

#	$Id: makefile,v 1.11 2013/02/03 Ning Exp $	

OBJS = AdjustStorage.o Aggregate.o AggregateRadiation.o	CalcAerodynamic.o \
CalcAvailableWater.o CalcBagnold.o CalcDistance.o CalcEffectiveKh.o \
CalcKhDry.o CalcKinViscosity.o CalcSafetyFactor.o CalcSatDensity.o CalcSnowAlbedo.o \
CalcSolar.o CalcTopoIndex.o \
CalcTotalWater.o CalcTransmissivity.o CalcWeights.o Calendar.o	     \
CanopyResistance.o ChannelState.o CheckOut.o CutBankGeometry.o	     \
DHSVMChannel.o Desorption.o DistSedDiams.o Draw.o EvalExponentIntegral.o \
EvapoTranspiration.o ExecDump.o FileIOBin.o FileIONetCDF.o Files.o   \
FinalMassBalance.o FindValue.o GetInit.o GetMetData.o InArea.o InitAggregated.o  \
InitArray.o InitConstants.o InitDump.o InitFileIO.o InitFineMaps.o   \
InitInterpolationWeights.o InitMetMaps.o InitMetSources.o	     \
InitModelState.o InitNetwork.o InitNewMonth.o InitParameters.o InitSedMap.o \
InitSedTables.o InitSnowMap.o InitTables.o InitTerrainMaps.o InitUnitHydrograph.o  \
InitXGraphics.o InterceptionStorage.o IsStationLocation.o LapseT.o LookupTable.o  \
MainDHSVM.o MainMWM.o MakeLocalMetData.o MassBalance.o MassEnergyBalance.o     \
MassRelease.o MaxRoadInfiltration.o NoEvap.o RadiationBalance.o	     \
ReadMetRecord.o ReadRadarMap.o ReportError.o ResetAggregate.o	     \
RootBrent.o Round.o RouteChannelSediment.o RouteRoad.o  RouteSubSurface.o   \
RouteSurface.o SatVaporPressure.o SensibleHeatFlux.o SeparateRadiation.o SizeOfNT.o \
SlopeAspect.o SnowInterception.o SnowMelt.o SnowPackEnergyBalance.o  \
SoilEvaporation.o StabilityCorrection.o StoreModelState.o	     \
SurfaceEnergyBalance.o UnsaturatedFlow.o VarID.o WaterTableDepth.o   \
channel.o channel_grid.o equal.o errorhandler.o globals.o tableio.o \
channel_complt.o CanopyShading.o

SRCS = $(OBJS:%.o=%.c)

HDRS = Calendar.h DHSVMChannel.h DHSVMerror.h brent.h channel.h	    \
channel_grid.h constants.h data.h errorhandler.h fifoNetCDF.h	    \
fifobin.h fileio.h functions.h getinit.h lookuptable.h massenergy.h \
rad.h settings.h sizeofnt.h slopeaspect.h snow.h soilmoisture.h	    \
tableio.h varid.h

OTHER = makefile tableio.lex

RCSDIR= RCS/
GET= co
REL=

 
DEFS =  -DHAVE_X11 -DHAVE_NETCDF
#possible DEFS -DHAVE_NETCDF -DHAVE_X11 -DSHOW_MET_ONLY -DSNOW_ONLY
CFLAGS =  -g -I/usr/X11R6/include -Wall  -I/usr/local/include/  $(DEFS) 

CC = cc
FLEX = /usr/bin/flex
LIBS = -lm -L/usr/X11R6/lib -lX11 -L/sw/lib -L/usr/local/lib -lnetcdf

# possible libs:   
#LIBS = -lm -L/usr/X11R6/lib -lX11 -L/sw/lib -L/usr/local/lib -lnetcdf

DHSVM: $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) -o DHSVM3.1.2 $(LIBS)

clean::
	rm -f DHSVM

library: libBinIO.a

BINIOOBJ = \
FileIOBin.o Files.o InitArray.o SizeOfNT.o Calendar.o \
ReportError.o

BINIOLIBOBJ = $(BINIOOBJ:%.o=libBinIO.a(%.o))

libBinIO.a: $(BINIOLIBOBJ)
	-ranlib $@

clean::
	rm -f libBinIO.a

# -------------------------------------------------------------
# rules for individual objects (created with make depend)
# -------------------------------------------------------------
AdjustStorage.o: AdjustStorage.c settings.h soilmoisture.h
Aggregate.o: Aggregate.c settings.h data.h Calendar.h DHSVMerror.h \
 functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
 constants.h
AggregateRadiation.o: AggregateRadiation.c settings.h data.h \
 Calendar.h massenergy.h
CalcAerodynamic.o: CalcAerodynamic.c DHSVMerror.h settings.h \
 constants.h functions.h data.h Calendar.h DHSVMChannel.h getinit.h \
 channel.h channel_grid.h
CalcAvailableWater.o: CalcAvailableWater.c settings.h soilmoisture.h
CalcBagnold.o: CalcBagnold.c DHSVMerror.h settings.h constants.h data.h \
 Calendar.h DHSVMChannel.h getinit.h channel.h channel_grid.h
CalcDistance.o: CalcDistance.c settings.h data.h Calendar.h \
 functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h
CalcEffectiveKh.o: CalcEffectiveKh.c settings.h constants.h \
 DHSVMerror.h functions.h data.h Calendar.h DHSVMChannel.h getinit.h \
 channel.h channel_grid.h
CalcKhDry.o: CalcKhDry.c settings.h functions.h data.h Calendar.h \
 DHSVMChannel.h getinit.h channel.h channel_grid.h
CalcKinViscosity.o: CalcKinViscosity.c functions.h data.h channel.h \
 DHSVMChannel.h settings.h Calendar.h getinit.h channel.h \
 channel_grid.h
CalcSafetyFactor.o: CalcSafetyFactor.c DHSVMerror.h settings.h \
 constants.h data.h Calendar.h
CalcSatDensity.o: CalcSatDensity.c settings.h functions.h constants.h
CalcSnowAlbedo.o: CalcSnowAlbedo.c settings.h constants.h data.h \
 Calendar.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h
CalcSolar.o: CalcSolar.c constants.h settings.h Calendar.h functions.h \
 data.h DHSVMChannel.h getinit.h channel.h channel_grid.h rad.h
CalcTopoIndex.o: CalcTopoIndex.c DHSVMerror.h constants.h data.h \
settings.h slopeaspect.h Calendar.h
CalcTotalWater.o: CalcTotalWater.c settings.h soilmoisture.h
CalcTransmissivity.o: CalcTransmissivity.c settings.h functions.h \
 data.h Calendar.h DHSVMChannel.h getinit.h channel.h channel_grid.h
CalcWeights.o: CalcWeights.c constants.h settings.h data.h Calendar.h \
 DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h
Calendar.o: Calendar.c settings.h functions.h data.h Calendar.h \
 DHSVMChannel.h getinit.h channel.h channel_grid.h DHSVMerror.h
CanopyShading.o: CanopyShading.c channel.h channel_grid.h constants.h \
 functions.h settings.h
CanopyResistance.o: CanopyResistance.c settings.h massenergy.h data.h \
 Calendar.h constants.h
channel_complt.o: channel_complt.c  functions.h errorhandler.h constants.h \
tableio.h settings.h
ChannelState.o: ChannelState.c settings.h data.h Calendar.h \
 DHSVMerror.h fileio.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h sizeofnt.h
CheckOut.o: CheckOut.c DHSVMerror.h settings.h data.h Calendar.h \
 functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
 constants.h
CutBankGeometry.o: CutBankGeometry.c settings.h soilmoisture.h
DHSVMChannel.o: DHSVMChannel.c constants.h getinit.h DHSVMChannel.h \
 settings.h data.h Calendar.h channel.h channel_grid.h DHSVMerror.h \
 functions.h errorhandler.h fileio.h
Desorption.o: Desorption.c settings.h massenergy.h data.h Calendar.h \
 constants.h
DistSedDiams.o: DistSedDiams.c data.h settings.h Calendar.h channel.h constants.h 
Draw.o: Draw.c settings.h data.h Calendar.h functions.h DHSVMChannel.h \
 getinit.h channel.h channel_grid.h snow.h
EvalExponentIntegral.o: EvalExponentIntegral.c settings.h data.h \
 Calendar.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h
EvapoTranspiration.o: EvapoTranspiration.c settings.h data.h \
 Calendar.h DHSVMerror.h massenergy.h constants.h
ExecDump.o: ExecDump.c settings.h data.h Calendar.h fileio.h \
 sizeofnt.h DHSVMerror.h functions.h DHSVMChannel.h getinit.h \
 channel.h channel_grid.h constants.h
FileIOBin.o: FileIOBin.c fifobin.h fileio.h sizeofnt.h settings.h \
 DHSVMerror.h
FileIONetCDF.o: FileIONetCDF.c
Files.o: Files.c settings.h data.h Calendar.h DHSVMerror.h functions.h \
 DHSVMChannel.h getinit.h channel.h channel_grid.h constants.h \
 fileio.h
FinalMassBalance.o: FinalMassBalance.c settings.h data.h Calendar.h \
 DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h
FindValue.o: FindValue.c data.h settings.h Calendar.h constants.h
GetInit.o: GetInit.c DHSVMerror.h fileio.h getinit.h
GetMetData.o: GetMetData.c settings.h data.h Calendar.h DHSVMerror.h \
 functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
 constants.h rad.h
InArea.o: InArea.c constants.h settings.h data.h Calendar.h
InitAggregated.o: InitAggregated.c settings.h data.h Calendar.h \
 DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h
InitArray.o: InitArray.c functions.h data.h settings.h Calendar.h \
 DHSVMChannel.h getinit.h channel.h channel_grid.h
InitConstants.o: InitConstants.c settings.h data.h Calendar.h fileio.h \
 DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h rad.h
InitDump.o: InitDump.c settings.h data.h Calendar.h DHSVMerror.h \
 fileio.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h sizeofnt.h varid.h
InitFileIO.o: InitFileIO.c fileio.h fifobin.h fifoNetCDF.h \
 DHSVMerror.h
InitFineMaps.o: InitFineMaps.c settings.h data.h Calendar.h \
 DHSVMerror.h fileio.h constants.h getinit.h varid.h sizeofnt.h
InitInterpolationWeights.o: InitInterpolationWeights.c settings.h \
 data.h Calendar.h DHSVMerror.h functions.h DHSVMChannel.h getinit.h \
 channel.h channel_grid.h constants.h
InitMetMaps.o: InitMetMaps.c settings.h constants.h data.h Calendar.h \
 DHSVMerror.h fileio.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h rad.h sizeofnt.h
InitMetSources.o: InitMetSources.c settings.h data.h Calendar.h \
 fileio.h DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h rad.h
InitModelState.o: InitModelState.c settings.h data.h Calendar.h \
 DHSVMerror.h fileio.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h sizeofnt.h soilmoisture.h varid.h
InitNetwork.o: InitNetwork.c constants.h settings.h data.h Calendar.h \
 DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h soilmoisture.h
InitNewMonth.o: InitNewMonth.c settings.h data.h Calendar.h \
 DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h fifobin.h fileio.h rad.h slopeaspect.h \
 sizeofnt.h
InitParameters.o: InitParameters.c settings.h data.h Calendar.h \
 fileio.h DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h rad.h
InitSedMap.o: InitSedMap.c data.h DHSVMerror.h 
InitSedTables.o: InitSedTables.c settings.h DHSVMerror.h Calendar.h \
 data.h constants.h fileio.h getinit.h
InitSnowMap.o: InitSnowMap.c settings.h data.h Calendar.h DHSVMerror.h \
 functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
 constants.h
InitTables.o: InitTables.c settings.h data.h Calendar.h DHSVMerror.h \
 functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
 constants.h fileio.h
InitTerrainMaps.o: InitTerrainMaps.c settings.h data.h Calendar.h \
 DHSVMerror.h fileio.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h sizeofnt.h slopeaspect.h varid.h
InitUnitHydrograph.o: InitUnitHydrograph.c settings.h constants.h \
 data.h Calendar.h DHSVMerror.h functions.h DHSVMChannel.h getinit.h \
 channel.h channel_grid.h fileio.h sizeofnt.h varid.h
InitXGraphics.o: InitXGraphics.c settings.h data.h Calendar.h \
 DHSVMerror.h
InterceptionStorage.o: InterceptionStorage.c settings.h data.h \
 Calendar.h DHSVMerror.h massenergy.h constants.h
IsStationLocation.o: IsStationLocation.c settings.h data.h Calendar.h \
 functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h
LapseT.o: LapseT.c settings.h data.h Calendar.h functions.h \
 DHSVMChannel.h getinit.h channel.h channel_grid.h
LookupTable.o: LookupTable.c lookuptable.h DHSVMerror.h
MainDHSVM.o: MainDHSVM.c settings.h constants.h data.h Calendar.h \
 DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h fileio.h
MainMWM.o: MainMWM.c settings.h Calendar.h getinit.h DHSVMerror.h \
 data.h fileio.h constants.h DHSVMChannel.h channel.h channel_grid.h
MakeLocalMetData.o: MakeLocalMetData.c settings.h data.h Calendar.h \
 snow.h DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h rad.h
MassBalance.o: MassBalance.c settings.h data.h Calendar.h DHSVMerror.h \
 functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
 constants.h
MassEnergyBalance.o: MassEnergyBalance.c settings.h data.h Calendar.h \
 DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h massenergy.h snow.h constants.h soilmoisture.h
MassRelease.o: MassRelease.c constants.h settings.h massenergy.h \
 data.h Calendar.h snow.h
MaxRoadInfiltration.o: MaxRoadInfiltration.c settings.h data.h \
 Calendar.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
 functions.h
NoEvap.o: NoEvap.c settings.h data.h Calendar.h massenergy.h
RadiationBalance.o: RadiationBalance.c settings.h data.h Calendar.h \
 DHSVMerror.h massenergy.h constants.h
ReadMetRecord.o: ReadMetRecord.c settings.h data.h Calendar.h \
 DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h
ReadRadarMap.o: ReadRadarMap.c settings.h data.h Calendar.h \
 DHSVMerror.h fileio.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h sizeofnt.h
ReportError.o: ReportError.c settings.h data.h Calendar.h DHSVMerror.h
ResetAggregate.o: ResetAggregate.c settings.h data.h Calendar.h \
 functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h
RootBrent.o: RootBrent.c settings.h brent.h massenergy.h data.h \
 Calendar.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h DHSVMerror.h
Round.o: Round.c functions.h data.h settings.h Calendar.h \
 DHSVMChannel.h getinit.h channel.h channel_grid.h DHSVMerror.h
RouteChannelSediment.o: data.h settings.h Calendar.h functions.h channel.h \
 constants.h DHSVMChannel.h getinit.h channel_grid.h DHSVMerror.h
RouteRoad.o: RouteRoad.c data.h settings.h Calendar.h DHSVMerror.h \
 functions.h channel.h DHSVMChannel.h constants.h channel_grid.h
RouteSubSurface.o: RouteSubSurface.c settings.h data.h Calendar.h \
 DHSVMerror.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h soilmoisture.h slopeaspect.h
RouteSurface.o: RouteSurface.c settings.h data.h Calendar.h \
 slopeaspect.h DHSVMerror.h functions.h DHSVMChannel.h getinit.h \
 channel.h channel_grid.h constants.h
SatVaporPressure.o: SatVaporPressure.c lookuptable.h
SensibleHeatFlux.o: SensibleHeatFlux.c settings.h data.h Calendar.h \
 DHSVMerror.h massenergy.h constants.h brent.h functions.h \
 DHSVMChannel.h getinit.h channel.h channel_grid.h
SeparateRadiation.o: SeparateRadiation.c settings.h rad.h
SizeOfNT.o: SizeOfNT.c DHSVMerror.h sizeofnt.h
SlopeAspect.o: SlopeAspect.c constants.h settings.h data.h Calendar.h \
 functions.h DHSVMChannel.h getinit.h channel.h channel_grid.h \
 slopeaspect.h DHSVMerror.h
SnowInterception.o: SnowInterception.c brent.h constants.h settings.h \
 massenergy.h data.h Calendar.h snow.h functions.h DHSVMChannel.h \
 getinit.h channel.h channel_grid.h
SnowMelt.o: SnowMelt.c brent.h constants.h settings.h massenergy.h \
 data.h Calendar.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h snow.h
SnowPackEnergyBalance.o: SnowPackEnergyBalance.c settings.h \
 constants.h massenergy.h data.h Calendar.h snow.h functions.h \
 DHSVMChannel.h getinit.h channel.h channel_grid.h
SoilEvaporation.o: SoilEvaporation.c settings.h DHSVMerror.h \
 massenergy.h data.h Calendar.h constants.h
StabilityCorrection.o: StabilityCorrection.c settings.h massenergy.h \
 data.h Calendar.h constants.h
StoreModelState.o: StoreModelState.c settings.h data.h Calendar.h \
 DHSVMerror.h fileio.h functions.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h constants.h sizeofnt.h varid.h
SurfaceEnergyBalance.o: SurfaceEnergyBalance.c settings.h massenergy.h \
 data.h Calendar.h constants.h
UnsaturatedFlow.o: UnsaturatedFlow.c constants.h settings.h \
 functions.h data.h Calendar.h DHSVMChannel.h getinit.h channel.h \
 channel_grid.h soilmoisture.h
VarID.o: VarID.c settings.h data.h Calendar.h DHSVMerror.h sizeofnt.h \
 varid.h
WaterTableDepth.o: WaterTableDepth.c settings.h soilmoisture.h
channel.o: channel.c errorhandler.h channel.h tableio.h settings.h
channel_grid.o: channel_grid.c channel_grid.h channel.h settings.h \
 data.h Calendar.h tableio.h errorhandler.h DHSVMChannel.h getinit.h
equal.o: equal.c functions.h data.h settings.h Calendar.h \
 DHSVMChannel.h getinit.h channel.h channel_grid.h
errorhandler.o: errorhandler.c errorhandler.h
globals.o: globals.c
tableio.o: tableio.c tableio.h errorhandler.h settings.h

tableio.c: tableio.lex
	$(FLEX) -Ptable_yy -o$@ $<

#clean::
#	rm -f tableio.c

# -------------------------------------------------------------
# sources
# -------------------------------------------------------------
sources: $(SRCS) $(HDRS) $(OTHER)

clean::
	rm -f $(OBJS)
	rm -f *~

# -------------------------------------------------------------
# tags 
# so we can find our way around
# -------------------------------------------------------------
tags: TAGS
TAGS: $(SRCS) $(HDRS)
	etags $(SRCS) $(HDRS)

clean::
	\rm -f TAGS


# -------------------------------------------------------------
# depend
# -------------------------------------------------------------
depend: .depend
.depend: $(SRCS)
	$(CC) $(CFLAGS) -MM $(SRCS) > $@

clean:: 
	rm -f .depend

