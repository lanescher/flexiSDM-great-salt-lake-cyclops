

# Install Rtools, if you do not have it already. Select the correct version
# of Rtools from this page: https://cran.r-project.org/bin/windows/Rtools/ 
# based on your version of R.
# to figure out what version of R you are using:
getRversion()

# Now restart R.


# Install the SpFut.flexiSDM package (this tutorial uses v1.1.2). Using it with
# a different version may cause errors.

# This package is not on CRAN, so you can't install it the way you would normally 
# install a package.

remotes::install_git(
  "https://code.usgs.gov/eastern-ecological-science-center/nearmi/SpFut-flexiSDM.git",
  ref = "1.1.2"
)
# When prompted, enter `1` to update all required packages
