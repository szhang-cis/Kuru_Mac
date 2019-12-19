import select
import socket
import os


class InfoServer(threading.Thread):
    def __init__(self):
        super(InfoServer, self).__init__()
        if not os.path.exists(os.path.expanduser(os.path.join("~", ".vulcanSock"))):
                os.mkdir(os.path.expanduser(os.path.join("~", ".vulcanSock")))
        self.server_address = os.path.expanduser(os.path.join("~", ".vulcanSock", "vulcan-" + str(pid) + ".sock"))
        try:
            os.unlink(self.server_address)
        except OSError:
            if os.path.exists(self.server_address):
                raise
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.alive = True
        self.pid = pid
        self.max_iter = 0
        self.iter_ = 0
        self.max_interval = 0
        self.interval = 0
        self.case_name = ""
        self.status = 0
        self.launch_time = 0
        self.elapsed = 0
        self.last_step_time = 0
        self.data_ready = False
    def setPid(self, pid):
        self.pid = pid
    def setCaseName(self, name):
        self.case_name = name
    def setStartTime(self, time):
        self.launch_time = time
    def run(self):
        self.sock.bind(self.server_address)
        self.sock.listen(5)
        while self.alive:
            (conn, address) = self.sock.accept()
            string = self.build_send_string()
            l = 0
            while l < len(string):
                sent = conn.send(string[l:])
                if sent == 0:
                    raise RuntimeError("Connection broken")
                l += sent
            while True:
                readable, writable, exceptional = select.select((conn,),(conn,),(conn,))
                for s in readable:
                    data = s.recv(10)
                    if data:
                        if data[0] == "K":
                            raise
                for s in writable:
                    if self.data_ready:
                        self.data_ready = False
                        l = 0
                        while l < len(self.outbuffer):
                            sent = s.send(self.outbuffer[l:])
                            if sent == 0:
                                raise RuntimeError("Connection broken")
                            l += sent
                for s  in exceptional:
                    pass        
    def stop(self):
        self.alive = False
    def baseInfo(self, max_iter, max_interval, case_name, launch_date):
        self.case_name = case_name
        self.max_iter = max_iter
        self.max_interval = max_interval
        self.launch_date = launch_date
    def newInfo(self, itera=False, interval=False, status=False, last_step_time=False):
        if last_step_time is not False:
            self.last_step_time = last_step_time
        if itera is not False:
            self.itera = itera
        if interval is not False:
            self.interval = interval
        if status is not False:
            self.status = status
        self.outbuffer = self.build_update_string()
        self.data_ready = True
    def build_send_string(self):
        #Format NP(PID 6 digit)I(maxiter 4 digit)N(maxinterval 3 digit)D(datetime Y-M-D H:M:S)C(casename 50 characters)
        return 'NP{pid:06}I{maxiter:04}N{maxinterval:03}D{datetime:%Y-%m-%d %H:%M:%S}C{casename:<50}'.format(pid=self.pid,maxiter=self.max_iter,maxinterval=self.max_interval,datetime=self.launch_date,casename=self.case_name)
    def build_update_string(self):
        # Format UI(iternum 4 digit)N(intervalnum 3 digit)S(s in 0 running, 1 error, 2 diverged, 3 converged)E(last step time in seconds 7 digits)
        string = 'UI{iternum:04}N{intervalnum:03}S{status}E{laststime:07}'.format(iternum = self.itera,intervalnum=self.interval,status=self.status,laststime=self.last_step_time)
        return string
    def address(self):
        return self.server_address
