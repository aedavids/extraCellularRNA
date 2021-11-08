
# capture output file
sink("hello.R.aedwip.out")
cat( "hello world \n")
print("good by")



# turn output capture off
sink()

# exit session do not save
q(save="no")
