

from mevis import *

def updateThreshold():
  """Updates threshold using the defined formula."""
  print "Calculating threshold"
  avnoise = ctx.field("StdDevNoise.totalMean").value
  modelrico = ctx.field("modelRico").value
  stdnoise = ctx.field("StdDevNoise.totalStdDev").value
  modelslope = ctx.field("modelSlope").value
  modeloffset = ctx.field("modelOffset").value
  tmin = (avnoise * modelrico) + (stdnoise * modelslope) + modeloffset 
  ctx.field("threshMin").value = tmin
  
#//# MeVis signature v1
#//# key: MFowDQYJKoZIhvcNAQEBBQADSQAwRgJBANEfsmYse2e1dRhkQ9AQbreCq9uxwzWLoGom13MNYmyfwoJqQOEXljLFAgw2eEjaT12G4CdqKWhRxh9ANP6n7GMCARE=:VI/mB8bT4u+mRtf/ru8yUQi8BzpaS3UeL2x62YxsUYnVqCWuLrVNLiukIIjnJMKQXlc8ezmgOIcVAV7pgvgKpQ==
#//# owner: MeVis
#//# date: 2013-06-24T21:32:11
#//# hash: bLEHIFu5kq+bqZpTtOVN0y3w9dJPhY1UZYZ37+dO638smoDvc3ftnwqDEo/j3QJp3WXZmhP64CUBSSCXwmD7Ug==
#//# MeVis end
