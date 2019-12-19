import curses
import time
import threading
import psutil
import datetime

class InfoServer(threading.Thread):
    def __init__(self):
        super(InfoServer, self).__init__()
        self.running = True
    def setPid(self, pid):
        self.pid = pid
        self.process = psutil.Process(self.pid)
    def setCaseName(self, name):
        self.case_name = name
    def setStartTime(self, time):
        self.launch_time = time
    def run(self):
        self.stdscr = curses.initscr()
        curses.start_color()
        curses.noecho()
        curses.cbreak()
        curses.curs_set(0)
        self.stdscr.keypad(1)
        curses.init_pair(1, curses.COLOR_YELLOW, curses.COLOR_BLACK)
        curses.init_pair(2, curses.COLOR_GREEN, curses.COLOR_GREEN)
        curses.init_pair(3, curses.COLOR_RED, curses.COLOR_WHITE)
        curses.init_pair(4, curses.COLOR_WHITE, curses.COLOR_BLACK)
        curses.init_pair(5, curses.COLOR_WHITE, curses.COLOR_YELLOW)
        y,x = self.stdscr.getmaxyx()
        self.stdscr.refresh()
        self.box(self.stdscr,0,0,50,13,1)
        self.stdscr.addch(9,0,curses.ACS_LTEE,curses.color_pair(1))
        self.stdscr.hline(9,1,curses.ACS_HLINE,29,curses.color_pair(1))
        self.stdscr.vline(1,30,curses.ACS_VLINE,12,curses.color_pair(1))
        self.stdscr.addch(9,30,curses.ACS_RTEE,curses.color_pair(1))
        self.stdscr.addch(0,30,curses.ACS_TTEE,curses.color_pair(1))
        self.stdscr.addch(13,30,curses.ACS_BTEE,curses.color_pair(1))
        self.drawWindow(self.stdscr, 50)
        self.button(self.stdscr, "Kill Process", 10,6, 5)
        self.Vbargraph(self.stdscr,1,32,5,10,00,1)
        self.stdscr.addstr(12,33,"Mem")
        self.stdscr.addstr(12,43,"CPU")
        while (self.running):
            self.stdscr.addstr(5,15,str(datetime.datetime.now()-self.launch_time).split(".")[0],curses.color_pair(4))
            self.stdscr.clrtoeol()
            self.Vbargraph(self.stdscr,1,32,5,10,self.process.memory_percent(),1)
            self.Vbargraph(self.stdscr,1,42,5,10,self.process.cpu_percent(),1)
            self.stdscr.refresh()
            
            time.sleep(0.1)
    def update(self,data):
        pass
    def stop(self):
        self.running = False
        curses.flash()
        curses.beep()
        curses.nocbreak()
        self.stdscr.keypad(0)
        curses.echo()
        curses.endwin()

    def Vbargraph(self,win,y,x,width,height,percent, color):
        bh = height - 2
        pbh = int(round(bh*(1-percent/100.0)))
        ph = int(round(bh*percent/100.0))
        self.box(win, y, x, width, height, color)
        if ph > 0:
            for i in range(1,width):
                win.vline(y+1+pbh,x+i," ",ph+1,curses.color_pair(2))
        win.addstr(y+height,x+1,"%3i%%" % percent)

    def button(self,win, text, y, x, color):
        width = len(text) + 5
        height = 2
        self.box(win, y, x, width, height, color)
        win.addstr(y+1,x+1,"  "+text+"  ",curses.color_pair(5))

    def box(self,win, y, x, width, height, color):
        win.hline(y,x+1,curses.ACS_HLINE,width-1,curses.color_pair(color))
        win.vline(y+1,x,curses.ACS_VLINE,height-1,curses.color_pair(color))
        win.vline(y+1,x+width,curses.ACS_VLINE,height-1,curses.color_pair(color))
        win.hline(y+height,x+1,curses.ACS_HLINE,width-1,curses.color_pair(color))
        win.addch(y,x, curses.ACS_ULCORNER, curses.color_pair(color))
        win.addch(y,x+width, curses.ACS_URCORNER, curses.color_pair(color))
        win.addch(y+height,x, curses.ACS_LLCORNER, curses.color_pair(color))
        win.addch(y+height,x+width, curses.ACS_LRCORNER, curses.color_pair(color))

    def drawWindow(self, win, x):
        text = "Vulcan"
        win.addch(0,((x-len(text))/2)-1, curses.ACS_RTEE, curses.color_pair(1))
        win.addch(0,((x+len(text))/2), curses.ACS_LTEE, curses.color_pair(1))
        win.addstr(0,(x-len(text))/2, text, curses.color_pair(5))
        win.addstr(1,1,"Case name: " + self.case_name, curses.color_pair(4))
        win.addstr(2,1,"Status: asdf", curses.color_pair(4))
        win.addstr(3,1,"PID: " + str(self.pid), curses.color_pair(4))
        win.addstr(4,1,"Start date: " + datetime.datetime.fromtimestamp(self.process.create_time()).strftime("%d-%m-%Y %H:%M:%S"), curses.color_pair(4))
        win.addstr(5,1,"Elapsed time: ", curses.color_pair(4))
        win.addstr(6,1,"Interval: 1/1", curses.color_pair(4))
        win.addstr(7,1,"Step: 1/1", curses.color_pair(4))
        win.addstr(8,1,"Last step time: asdf", curses.color_pair(4))

def main():    
    stdscr = curses.initscr()
    curses.start_color()
    curses.noecho()
    curses.cbreak()
    curses.curs_set(0)
    stdscr.keypad(1)
    curses.init_pair(1, curses.COLOR_YELLOW, curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_GREEN, curses.COLOR_GREEN)
    curses.init_pair(3, curses.COLOR_RED, curses.COLOR_WHITE)
    curses.init_pair(4, curses.COLOR_WHITE, curses.COLOR_BLACK)
    curses.init_pair(5, curses.COLOR_WHITE, curses.COLOR_YELLOW)
    y,x = stdscr.getmaxyx()
    stdscr.refresh()
    box(stdscr,0,0,50,13,1)
    stdscr.addch(9,0,curses.ACS_LTEE,curses.color_pair(1))
    stdscr.hline(9,1,curses.ACS_HLINE,29,curses.color_pair(1))
    stdscr.vline(1,30,curses.ACS_VLINE,12,curses.color_pair(1))
    stdscr.addch(9,30,curses.ACS_RTEE,curses.color_pair(1))
    stdscr.addch(0,30,curses.ACS_TTEE,curses.color_pair(1))
    stdscr.addch(13,30,curses.ACS_BTEE,curses.color_pair(1))
    drawWindow(stdscr, 50)
    button(stdscr, "Kill Process", 10,6, 5)
    Vbargraph(stdscr,1,32,5,10,00,1)
    stdscr.addstr(12,33,"Mem")
    stdscr.addstr(12,43,"CPU")
    for i in range(0,100,1):
        Vbargraph(stdscr,1,32,5,10,i,1)
        Vbargraph(stdscr,1,42,5,10,i,1)
        stdscr.refresh()
        time.sleep(0.1)

    curses.flash()
    curses.beep()
    curses.nocbreak()
    stdscr.keypad(0)
    curses.echo()
    curses.endwin()

if __name__ == "__main__":
    main()