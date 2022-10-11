#!/bin/csh -ef

# blatJob.csh - script to run the blat operation during the
#           - same species liftOver procedure
#           - ASSUMING current working directory is the ./run.blat/
#           - work directory so that targetList and queryListIn references
#           - will be correct at ./tParts/partXXX.lst and ./qParts/partXXX.lst
#           - and 11.ooc exists here for the target genome sequence

set targetList = $1
set queryListIn = $2
set outPsl = $3

if ($targetList:e == "lst") set targetList = $targetList
if ($queryListIn:e == "lst") set queryListIn = $queryListIn

# Use local disk for output, and move the final result to $outPsl
# when done, to minimize I/O.
set tmpDir = `mktemp -d -p /dev/shm doSame.blat.XXXXXX`
pushd $tmpDir

# We might get a .lst or a 2bit spec here -- convert to (list of) 2bit spec:
if ($queryListIn:e == "lst") then
  set specList = `cat $queryListIn`
else
  set specList = $queryListIn
endif

# Further partition the query spec(s) into 5000 coord ranges, building up
# a .lst of 2bit specs for blat and a .lft liftUp spec for the results:
cp /dev/null reSplitQuery.lst
cp /dev/null query.lft
foreach spec ($specList)
  set file  = `echo $spec | awk -F: '{print $1;}'`
  set seq   = `echo $spec | awk -F: '{print $2;}'`
  set range = `echo $spec | awk -F: '{print $3;}'`
  set start = `echo $range | awk -F- '{print $1;}'`
  set end   = `echo $range | awk -F- '{print $2;}'`
  if (! -e q.sizes) twoBitInfo $file q.sizes
  set seqSize = `awk '$1 == "'$seq'" {print $2;}' q.sizes`
  set chunkEnd = '0'
  while ($chunkEnd < $end)
    set chunkEnd = `echo $start 5000 | awk '{print $1+$2}'`
    if ($chunkEnd > $end) set chunkEnd = $end
    set chunkSize = `echo $chunkEnd $start | awk '{print $1-$2}'`
    echo $file\:$seq\:$start-$chunkEnd >> reSplitQuery.lst
    if (($start == 0) && ($chunkEnd == $seqSize)) then
      echo "$start      $seq    $seqSize        $seq    $seqSize" >> query.lft
    else
      echo "$start      $seq"":$start-$chunkEnd $chunkSize      $seq    $seqSize" >> query.lft
    endif
    set start = `echo $chunkEnd 500 | awk '{print $1-$2}'`
  end
end

# Align unsplit target sequence(s) to .lst of 2bit specs for 5000 chunks
# of query:
blat $targetList reSplitQuery.lst tmpUnlifted.psl \
  -tileSize=11 -ooc=/mnt/DATA4/kevhu/programs/liftOverSameSpecies/run.blat/11.ooc \
  -minScore=100 -minIdentity=98 -fastMap -noHead

# Lift query coords back up:
liftUp -pslQ -nohead tmpOut.psl query.lft warn tmpUnlifted.psl

# Move final result into place:
mv tmpOut.psl $outPsl

popd
rm -rf $tmpDir
