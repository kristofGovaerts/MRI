def previousLabel():
  ctx.field(Calculator.i1).value = ctx.field(RunPythonScript.out1).value
  
def nextLabel():
  ctx.field(Calculator.i1).value = ctx.field(RunPythonScript.out2).value