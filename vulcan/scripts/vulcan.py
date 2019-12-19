#!/usr/bin/python
# -*- coding: utf-8 -*-
from multiprocessing import Process, cpu_count
import os
import shutil
import subprocess
import select
import datetime


class Vulcan(Process):
    def __init__(self, res_directory, Arguments):
        Process.__init__(self)
        self.pathToVulcan = "/usr/bin/vulcan-bin"
        self.arguments = Arguments
        self.res_directory = res_directory

#       Definicion modelo de ejecucion
        self.modelo = None  # Antiguo (orden rigido) o Nuevo (orden flexible)
        if (len(self.arguments[0]) > 1) or (len(self.arguments[0]) == 1 and len(self.arguments) == 2):  # Modelo Antiguo
            self.modelo = "Antiguo"
        elif len(self.arguments) % 2 == 0:
            self.modelo = "Nuevo"
        else:
            raise
        self.filename = {
            't': None,
            'f': None,
            'm': None,
            's': None,
            'c': None,
        }
        self.dat_directory = os.path.dirname(os.path.abspath(self.arguments[1]))
        scratch_name = os.path.splitext(os.path.basename(self.arguments[1]))[0]
        if os.path.splitext(self.arguments[1].split("-")[-1])[0] in ['t', 'f', 'm', 's', 'c']:
            self.res_directory = os.path.join(res_directory, "-".join(scratch_name.split("-")[:-1]))
        else:
            self.res_directory = os.path.join(res_directory, scratch_name)

        if self.modelo == "Antiguo":
            vulcan_type = self.arguments[0]
            filename_list = [os.path.splitext(os.path.basename(x))[0] for x in self.arguments[1:]]  # Obtiene el nombre del archivo sin extension sin directorios
            if vulcan_type == "t":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-t.O2")
                self.filename['t'] = filename_list[0]
                self.res_directory += '-t'
            if vulcan_type == "f":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-f.O2")
                self.filename['f'] = filename_list[0]
                self.res_directory += '-f'
            if vulcan_type == "m":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-m.O2")
                self.filename['m'] = filename_list[0]
                self.res_directory += '-m'
            if vulcan_type == "s":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-s.O2")
                self.filename['s'] = filename_list[0]
                self.res_directory += '-s'
            if vulcan_type == "tm":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-tm.O2")
                self.filename['t'] = filename_list[0]
                self.filename['m'] = filename_list[1]
                self.filename['c'] = filename_list[2]
                self.res_directory += '-tm'
            if vulcan_type == "tf":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-tf.O2")
                self.filename['t'] = filename_list[0]
                self.filename['f'] = filename_list[1]
                self.res_directory += '-tf'
            if vulcan_type == "ts":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-ts.O2")
                self.filename['t'] = filename_list[0]
                self.filename['s'] = filename_list[1]
                self.res_directory += '-ts'
            if vulcan_type == "mf":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-mf.O2")
                self.filename['m'] = filename_list[0]
                self.filename['f'] = filename_list[1]
                self.res_directory += '-mf'
            if vulcan_type == "ms":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-ms.O2")
                self.filename['m'] = filename_list[0]
                self.filename['s'] = filename_list[1]
                self.res_directory += '-ms'
            if vulcan_type == "tms":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-tms.O2")
                self.filename['t'] = filename_list[0]
                self.filename['m'] = filename_list[1]
                self.filename['s'] = filename_list[2]
                self.filename['c'] = filename_list[3]
                self.res_directory += '-tms'
            if vulcan_type == "tfm":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-tfm.O2")
                self.filename['t'] = filename_list[0]
                self.filename['f'] = filename_list[1]
                self.filename['m'] = filename_list[2]
                self.res_directory += '-tfm'
            if vulcan_type == "tfms":
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-tfms.O2")
                self.filename['t'] = filename_list[0]
                self.filename['f'] = filename_list[1]
                self.filename['m'] = filename_list[2]
                self.filename['s'] = filename_list[3]
                self.filename['c'] = filename_list[4]
                self.res_directory += '-tfms'
        elif self.modelo == "Nuevo":
            order = {
                "t": 1,
                "f": 2,
                "m": 3,
                "s": 4,
                "c": 5,
            }
            rundict = {self.arguments[i]: os.path.splitext(os.path.basename(self.arguments[i+1]))[0] for i in range(0,len(self.arguments),2)}
            for key, value in rundict.iteritems():
                self.filename[key] = value
            sortedlist = sorted(rundict.keys(), key=lambda x: order[x[0]])
            if (
                (['t', 'm'] == sortedlist[:-1] and 'c' not in sortedlist) or
                (['t', 'm', 's'] == sortedlist[:-1] and 'c' not in sortedlist) or
                (['t', 'f', 'm', 's'] == sortedlist[:-1] and 'c' not in sortedlist)
            ):
                raise NameError("Missing coupling file")
            if 'c' not in sortedlist:
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-" + ''.join(sortedlist) + ".O2")
                self.res_directory += '-' + ''.join(sortedlist)
            else:
                self.vulcan_exec = os.path.join(self.pathToVulcan, "Vulcan-" + ''.join(sortedlist[:-1]) + ".O2")
                self.res_directory += '-' + ''.join(sortedlist[:-1])

        for key, item in self.filename.iteritems():
            if item is not None:
                if not os.path.isfile(os.path.join(self.dat_directory, item + ".dat")):
                    raise NameError(os.path.join(self.dat_directory, item + ".dat"))

#       Variables de Entorno
        self.env = {}
        self.env_temp = {}
        self.env_res = {}

        self.general_name = os.path.basename(self.res_directory)

        self.env["GFORTRAN_UNBUFFERED_ALL"] = "1"
        self.env['LD_LIBRARY_PATH'] = '/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64:/opt/intel/composer_xe_2013_sp1.3.174/mpirt/lib/intel64:/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64:/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/intel64'

#       General
        self.env_res["FOR500"] = os.path.join(self.res_directory, self.general_name + ".tim")  # false unit representing stdout
        self.env_res["FOR501"] = os.path.join(self.res_directory, self.general_name + ".err")  # false unit representing stderr
        self.env_temp["FOR502"] = os.path.join(self.res_directory, self.general_name + ".inf")  # Status output THIS IS A FIFO!!

#       Vulcan-M
        if self.filename['m'] is not None:
#           Dat directory
            self.env["FOR105"] = os.path.join(self.dat_directory, self.filename['m'] + ".dat")
            self.env["FOR140"] = os.path.join(self.dat_directory, self.filename['m'] + ".geo")
            self.env["FOR141"] = os.path.join(self.dat_directory, self.filename['m'] + ".set")
            self.env["FOR142"] = os.path.join(self.dat_directory, self.filename['m'] + ".mat")
            self.env["FOR143"] = os.path.join(self.dat_directory, self.filename['m'] + ".ini")
            self.env["FOR144"] = os.path.join(self.dat_directory, self.filename['m'] + ".loa")
            self.env["FOR145"] = os.path.join(self.dat_directory, self.filename['m'] + ".fix")
            self.env["FOR146"] = os.path.join(self.dat_directory, self.filename['m'] + ".ini1")
            self.env["FOR147"] = os.path.join(self.dat_directory, self.filename['m'] + ".tun")
            self.env["FOR148"] = os.path.join(self.dat_directory, self.filename['m'] + ".con")
            self.env["FOR149"] = os.path.join(self.dat_directory, self.filename['m'] + ".act")
#           Scratch directory
            self.env_res["FOR101"] = os.path.join(self.res_directory, self.filename['m'] + ".dts")
            self.env_res["FOR102"] = os.path.join(self.res_directory, self.filename['m'] + ".sol")
            self.env_res["FOR103"] = os.path.join(self.res_directory, self.filename['m'] + ".fro")
            self.env_res["FOR104"] = os.path.join(self.res_directory, self.filename['m'] + ".frhk")
            self.env_res["FOR106"] = os.path.join(self.res_directory, self.filename['m'] + ".pri")
            self.env_res["FOR107"] = os.path.join(self.res_directory, self.filename['m'] + ".res")
            self.env_res["FOR108"] = os.path.join(self.res_directory, self.filename['m'] + ".sol2")
            self.env_res["FOR109"] = os.path.join(self.res_directory, self.filename['m'] + ".fro2")
            self.env_res["FOR110"] = os.path.join(self.res_directory, self.filename['m'] + ".pos")
            self.env_res["FOR111"] = os.path.join(self.res_directory, self.filename['m'] + ".rst")
            self.env_res["FOR112"] = os.path.join(self.res_directory, self.filename['m'] + ".bfgs")
            self.env_res["FOR113"] = os.path.join(self.res_directory, self.filename['m'] + ".pipe")
            self.env_res["FOR114"] = os.path.join(self.res_directory, self.filename['m'] + ".pan")
            self.env_res["FOR139"] = os.path.join(self.res_directory, self.filename['m'] + ".fan")
            self.env_res["FOR191"] = os.path.join(self.res_directory, self.filename['m'] + ".cur1")
            self.env_res["FOR192"] = os.path.join(self.res_directory, self.filename['m'] + ".cur2")
            self.env_res["FOR193"] = os.path.join(self.res_directory, self.filename['m'] + ".cur3")
            self.env_res["FOR194"] = os.path.join(self.res_directory, self.filename['m'] + ".cur4")
            self.env_res["FOR195"] = os.path.join(self.res_directory, self.filename['m'] + ".cur5")
            self.env_res["FOR196"] = os.path.join(self.res_directory, self.filename['m'] + ".cur6")
            self.env_res["FOR197"] = os.path.join(self.res_directory, self.filename['m'] + ".cur7")
            self.env_res["FOR198"] = os.path.join(self.res_directory, self.filename['m'] + ".cur8")
            self.env_res["FOR199"] = os.path.join(self.res_directory, self.filename['m'] + ".cur9")
            self.env_res["FOR200"] = os.path.join(self.res_directory, self.filename['m'] + ".cur10")

#       Vulcan-T
        if self.filename['t'] is not None:
#           Dat directory
            self.env["FOR225"] = os.path.join(self.dat_directory, self.filename['t'] + ".dat")
            self.env["FOR250"] = os.path.join(self.dat_directory, self.filename['t'] + ".geo")
            self.env["FOR251"] = os.path.join(self.dat_directory, self.filename['t'] + ".set")
            self.env["FOR252"] = os.path.join(self.dat_directory, self.filename['t'] + ".mat")
            self.env["FOR253"] = os.path.join(self.dat_directory, self.filename['t'] + ".ini")
            self.env["FOR254"] = os.path.join(self.dat_directory, self.filename['t'] + ".loa")
            self.env["FOR255"] = os.path.join(self.dat_directory, self.filename['t'] + ".fix")
            self.env["FOR256"] = os.path.join(self.dat_directory, self.filename['t'] + ".adv")
            self.env["FOR259"] = os.path.join(self.dat_directory, self.filename['t'] + ".act")
            self.env["FOR282"] = os.path.join(self.dat_directory, self.filename['t'] + ".str")
#           Scratch directory
            self.env_res["FOR221"] = os.path.join(self.res_directory, self.filename['t'] + ".dts")
            self.env_res["FOR222"] = os.path.join(self.res_directory, self.filename['t'] + ".sol")
            self.env_res["FOR223"] = os.path.join(self.res_directory, self.filename['t'] + ".fro")
            self.env_res["FOR224"] = os.path.join(self.res_directory, self.filename['t'] + ".frhk")
            self.env_res["FOR226"] = os.path.join(self.res_directory, self.filename['t'] + ".pri")
            self.env_res["FOR227"] = os.path.join(self.res_directory, self.filename['t'] + ".res")
            self.env_res["FOR228"] = os.path.join(self.res_directory, self.filename['t'] + ".sol2")
            self.env_res["FOR229"] = os.path.join(self.res_directory, self.filename['t'] + ".fro2")
            self.env_res["FOR230"] = os.path.join(self.res_directory, self.filename['t'] + ".pos")
            self.env_res["FOR231"] = os.path.join(self.res_directory, self.filename['t'] + ".rst")
            self.env_res["FOR232"] = os.path.join(self.res_directory, self.filename['t'] + ".bfgs")
            self.env_res["FOR233"] = os.path.join(self.res_directory, self.filename['t'] + ".pipe")
            self.env_res["FOR234"] = os.path.join(self.res_directory, self.filename['t'] + ".pan")
            self.env_res["FOR238"] = os.path.join(self.res_directory, self.filename['t'] + ".fan")
            self.env_res["FOR291"] = os.path.join(self.res_directory, self.filename['t'] + ".cur1")
            self.env_res["FOR292"] = os.path.join(self.res_directory, self.filename['t'] + ".cur2")
            self.env_res["FOR293"] = os.path.join(self.res_directory, self.filename['t'] + ".cur3")
            self.env_res["FOR294"] = os.path.join(self.res_directory, self.filename['t'] + ".cur4")
            self.env_res["FOR295"] = os.path.join(self.res_directory, self.filename['t'] + ".cur5")
            self.env_res["FOR296"] = os.path.join(self.res_directory, self.filename['t'] + ".cur6")
            self.env_res["FOR297"] = os.path.join(self.res_directory, self.filename['t'] + ".cur7")
            self.env_res["FOR298"] = os.path.join(self.res_directory, self.filename['t'] + ".cur8")
            self.env_res["FOR299"] = os.path.join(self.res_directory, self.filename['t'] + ".cur9")
            self.env_res["FOR300"] = os.path.join(self.res_directory, self.filename['t'] + ".cur10")

#       Vulcan-F
        if self.filename['f'] is not None:
#           Dat directory
            self.env["FOR092"] = os.path.join(self.dat_directory, self.filename['f'] + ".dat")
            self.env["FOR084"] = os.path.join(self.dat_directory, self.filename['f'] + ".str")
            self.env["FOR091"] = os.path.join(self.dat_directory, self.filename['f'] + ".fri")  # Frictional conditions file
            self.env["FOR080"] = os.path.join(self.dat_directory, self.filename['f'] + ".rhs")  # Subincrementation of RHS file
            self.env["FOR057"] = os.path.join(self.dat_directory, self.filename['f'] + ".vbo")  # Variable boundary contition file
            self.env["FOR058"] = os.path.join(self.dat_directory, self.filename['f'] + ".rvb")  # Restart variable boundary contition file
#           Scratch directory
            self.env_res["FOR093"] = os.path.join(self.res_directory, self.filename['f'] + ".inc")  # include file
            self.env_res["FOR094"] = os.path.join(self.res_directory, self.filename['f'] + ".log")  # include file
            self.env_res["FOR095"] = os.path.join(self.res_directory, self.filename['f'] + ".res")  # include file
            self.env_res["FOR096"] = os.path.join(self.res_directory, self.filename['f'] + ".pos")  # include file
            self.env_res["FOR097"] = os.path.join(self.res_directory, self.filename['f'] + ".rst")  # include file
            self.env_res["FOR098"] = os.path.join(self.res_directory, self.filename['f'] + ".rsf")  # include file
            self.env_res["FOR081"] = os.path.join(self.res_directory, self.filename['f'] + ".con")  # Contour position output file

            self.env_temp["FOR085"] = os.path.join(self.res_directory, self.filename['f'] + ".dts")
            self.env_temp["FOR086"] = os.path.join(self.res_directory, self.filename['f'] + ".sc1")  # solver file 1
            self.env_temp["FOR087"] = os.path.join(self.res_directory, self.filename['f'] + ".sc2")  # solver file 2
            self.env_temp["FOR088"] = os.path.join(self.res_directory, self.filename['f'] + ".sc3")  # solver file 3
            self.env_temp["FOR089"] = os.path.join(self.res_directory, self.filename['f'] + ".sc4")  # solver file 4
            self.env_temp["FOR099"] = os.path.join(self.res_directory, self.filename['f'] + ".sc5")  # solver file 5

#       Vulcan-S
        if self.filename['s'] is not None:
#           Dat directory
            self.env["FOR315"] = os.path.join(self.dat_directory, self.filename['s'] + ".dat")
            self.env["FOR350"] = os.path.join(self.dat_directory, self.filename['s'] + ".geo")
            self.env["FOR351"] = os.path.join(self.dat_directory, self.filename['s'] + ".set")
            self.env["FOR352"] = os.path.join(self.dat_directory, self.filename['s'] + ".mat")
            self.env["FOR353"] = os.path.join(self.dat_directory, self.filename['s'] + ".ini")
            self.env["FOR354"] = os.path.join(self.dat_directory, self.filename['s'] + ".loa")
            self.env["FOR355"] = os.path.join(self.dat_directory, self.filename['s'] + ".fix")
            self.env["FOR356"] = os.path.join(self.dat_directory, self.filename['s'] + ".adv")
            self.env["FOR359"] = os.path.join(self.dat_directory, self.filename['s'] + ".act")
            self.env["FOR382"] = os.path.join(self.dat_directory, self.filename['s'] + ".str")
#           Scratch directory
            self.env_res["FOR316"] = os.path.join(self.res_directory, self.filename['s'] + ".res")
            self.env_res["FOR317"] = os.path.join(self.res_directory, self.filename['s'] + ".sol")
            self.env_res["FOR318"] = os.path.join(self.res_directory, self.filename['s'] + ".frhk")
            self.env_res["FOR319"] = os.path.join(self.res_directory, self.filename['s'] + ".pos")
            self.env_res["FOR320"] = os.path.join(self.res_directory, self.filename['s'] + ".dts")
            self.env_res["FOR321"] = os.path.join(self.res_directory, self.filename['s'] + ".fro")
            self.env_res["FOR322"] = os.path.join(self.res_directory, self.filename['s'] + ".pri")
            self.env_res["FOR323"] = os.path.join(self.res_directory, self.filename['s'] + ".sol2")
            self.env_res["FOR324"] = os.path.join(self.res_directory, self.filename['s'] + ".fro2")
            self.env_res["FOR325"] = os.path.join(self.res_directory, self.filename['s'] + ".rst")
            self.env_res["FOR326"] = os.path.join(self.res_directory, self.filename['s'] + ".bfgs")
            self.env_res["FOR327"] = os.path.join(self.res_directory, self.filename['s'] + ".pipe")
            self.env_res["FOR328"] = os.path.join(self.res_directory, self.filename['s'] + ".pan")
            self.env_res["FOR329"] = os.path.join(self.res_directory, self.filename['s'] + ".fan")
            self.env_res["FOR391"] = os.path.join(self.res_directory, self.filename['s'] + ".cur1")
            self.env_res["FOR392"] = os.path.join(self.res_directory, self.filename['s'] + ".cur2")
            self.env_res["FOR393"] = os.path.join(self.res_directory, self.filename['s'] + ".cur3")
            self.env_res["FOR394"] = os.path.join(self.res_directory, self.filename['s'] + ".cur4")
            self.env_res["FOR395"] = os.path.join(self.res_directory, self.filename['s'] + ".cur5")
            self.env_res["FOR396"] = os.path.join(self.res_directory, self.filename['s'] + ".cur6")
            self.env_res["FOR397"] = os.path.join(self.res_directory, self.filename['s'] + ".cur7")
            self.env_res["FOR398"] = os.path.join(self.res_directory, self.filename['s'] + ".cur8")
            self.env_res["FOR399"] = os.path.join(self.res_directory, self.filename['s'] + ".cur9")
            self.env_res["FOR400"] = os.path.join(self.res_directory, self.filename['s'] + ".cur10")

#       Coupling units
        if self.filename['c'] is not None:
#           Dat directory
            self.env["FOR435"] = os.path.join(self.dat_directory, self.filename['c'] + ".dat")
#           Scratch directory
            self.env["FOR436"] = os.path.join(self.res_directory, self.filename['c'] + ".res")

        if os.path.isdir(self.res_directory):
            for key, value in self.env_res.iteritems():
                try:
                    os.remove(value)
                except OSError:
                    continue
        else:
            os.makedirs(self.res_directory)

        bak = ["dat", "geo", "set", "mat", "ini", "loa", "fix", "ini1", "tun", "con", "act", "adv", "str", "fri", "rhs", "vbo", "rvb"]
        for item in bak:
            try:
                for fname in [j for j in self.filename.itervalues() if j is not None]:
                    shutil.copy((os.path.join(self.dat_directory, fname + "." + item)), self.res_directory)
            except IOError:
                continue

    def run(self, InfoServer=None):
        self.env.update(self.env_temp)
        self.env.update(self.env_res)
        try:
            os.remove(self.env_temp["FOR502"])
        except OSError:
            pass
        os.mkfifo(self.env_temp["FOR502"])
        p = subprocess.Popen(self.vulcan_exec, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.res_directory, env=self.env)
        fifo_ = open(self.env_temp["FOR502"], "r", 0)
        fifo_fd = fifo_.fileno()
        print " "
        print "Ejecutando Vulcan Script PID:", os.getpid()
        print "Ejecutando Vulcan.O2     PID:", p.pid
        if InfoServer is not None:
            InfoServer.setPid(p.pid)
            InfoServer.setCaseName(self.general_name)
            InfoServer.setStartTime(datetime.datetime.now())
            InfoServer.start()
        buf = ""
        while p.poll() is None:
            rlist, wlist, xlist = select.select((fifo_fd,), (), (),100)
            for item in rlist:
                c = os.read(item,1)
                buf += c
                if c == "\n":
                    print buf
                    if InfoServer is not None:
                        InfoServer.update(buf[:-1])
                    buf = ""
        if InfoServer is not None:
            InfoServer.stop()
        print "Vulcan se terminó de ejecutar"
        self.stdout = p.stdout
        self.stderr = p.stderr
        fifo_.close()
        for key, item in self.env_temp.iteritems():
            try:
                os.remove(item)
            except OSError:
                continue
        self.returnValue = p.returncode
        tim = open(self.env['FOR500'], 'w')
        err = open(self.env['FOR501'], 'w')
        tim.writelines(self.stdout.readlines())
        err.writelines(self.stderr.readlines())
        tim.flush()
        tim.close()
        err.flush()
        err.close()
