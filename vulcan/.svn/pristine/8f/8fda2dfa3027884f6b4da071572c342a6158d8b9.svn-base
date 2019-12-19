#!/usr/bin/env python

import pysvn, subprocess, os

repository = "/home/vulcan-build/vulcan"
prerepository = "/home/vulcan-build"
version_major = 1
version_minor = 2
client = pysvn.Client()
rev = client.update(repository)[0]

env = os.environ.copy()
env["DEB_FFLAGS_MAINT_APPEND"] = "-fopenmp -cpp"
env["MKLROOT"] = "/opt/intel/composer_xe_2015.0.090/mkl"
env["LD_LIBRARY_PATH"] = "/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64:/opt/intel/composer_xe_2015.0.090/mpirt/lib/intel64:/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64:/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64"
env["NLSPATH"]= "/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64/locale/%l_%t/%N:/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64/locale/%l_%t/%N:/opt/intel/composer_xe_2015.0.090/debugger/gdb/intel64_mic/share/locale/%l_%t/%N:/opt/intel/composer_xe_2015.0.090/debugger/gdb/intel64/share/locale/%l_%t/%N"
env["MPM_LAUNCHER"]= "/opt/intel/composer_xe_2015.0.090/debugger/mpm/bin/start_mpm.sh"
env["MIC_LIBRARY_PATH"]= "/opt/intel/composer_xe_2015.0.090/compiler/lib/mic:/opt/intel/composer_xe_2015.0.090/mpirt/lib/mic"
env["MIC_LD_LIBRARY_PATH"]= "/opt/intel/composer_xe_2015.0.090/compiler/lib/mic:/opt/intel/composer_xe_2015.0.090/mpirt/lib/mic:/opt/intel/composer_xe_2015.0.090/compiler/lib/mic:/opt/intel/composer_xe_2015.0.090/mkl/lib/mic"
env["LIBRARY_PATH"]= "/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64:/opt/intel/composer_xe_2015.0.090/compiler/lib/intel64:/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64"
env["INTEL_LICENSE_FILE"]= "/opt/intel/composer_xe_2015.0.090/licenses:/opt/intel/licenses:/home/vulcan-build/intel/licenses"
env["INTEL_PYTHONHOME"]= "/opt/intel/composer_xe_2015.0.090/debugger/python/intel64/"
env["INCLUDE"]= "/opt/intel/composer_xe_2015.0.090/mkl/include"
env["INFOPATH"]= "/opt/intel/composer_xe_2015.0.090/debugger/gdb/intel64/share/info/:/opt/intel/composer_xe_2015.0.090/debugger/gdb/intel64_mic/share/info/"

print "Cleaning..."
p = subprocess.Popen(["debian/rules", "clean"], env=env, cwd=repository, shell=True)
p.wait()
print "Making orig tarball"
p = subprocess.Popen(["tar", "czf", "vulcan_"+str(version_major)+"."+str(version_minor)+".orig.tar.gz", "vulcan"], cwd=prerepository)
p.wait()

p = subprocess.Popen(["dch","-v",str(version_major)+"."+str(version_minor)+"-r"+str(rev.number),client.log(repository,limit=0)[0]['message']], cwd=repository, env=env)
p.wait()

p = subprocess.call(["dpkg-buildpackage", "-us", "-uc"], cwd=repository, env=env)
#p.wait()

for dist in ["precise", "trusty", "wheezy", "jessie"]:
	p = subprocess.call(["reprepro", "-b/srv/vulcan-repo", "includedeb", dist, "vulcan_"+str(version_major)+"."+str(version_minor)+"-r"+str(rev.number)+"_amd64.deb"], env=env, cwd=prerepository)
