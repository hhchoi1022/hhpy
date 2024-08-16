import serial
import time
import numpy as np

class AntSerial():
    def __init__(self, strPort = None, nBaudrate = 9600, strParity = 'N', nDatabits = 8, nStopbits = 1, nTimeout = 1, bXonxoff = False):

        self.com = serial.Serial()

        # serial communication parameter
        self.com.port = strPort
        self.com.baudrate = nBaudrate
        self.com.parity = strParity
        self.com.bytesize = nDatabits
        self.com.stopbits = nStopbits
        self.com.timeout = nTimeout
        self.com.xonxoff = bXonxoff

        # communication command/response rule
        self.sendSuffix = b'\r\n'
        self.receiveSuffix = b'\n'
        self.receiveMaxLength = 96
        self.waitingTime = 3.0        
        self.receiveFrontIndex = None
        self.receiveBackIndex = None

    def SetParameter(self, strPort = None, nBaudrate = 9600, strParity = 'N', nDatabits = 8, nStopbits = 1, nTimeout = 1):
        try:
            self.com.port = strPort
            self.com.baudrate = nBaudrate
            self.com.parity = strParity
            self.com.bytesize = nDatabits
            self.com.stopbits = nStopbits
            self.com.timeout = nTimeout
            return 1
        except:
            return 0

    def CommunicationRule(self, bCommandSuffix = b'\r\n', bResponseSuffix = b'\n', nResponseMaxLength = 96, fwaitingTime = 3.0,
    nResponseSliceFrontIndex = None, nResponseSliceBackIndex = None):        
        # print(bCommandSuffix, bResponseSuffix, nResponseMaxLength, fwaitingTime, nResponseSliceFrontIndex, nResponseSliceBackIndex)
        try:
            self.sendSuffix = bCommandSuffix
            self.receiveSuffix = bResponseSuffix

            self.receiveMaxLength = nResponseMaxLength
            self.waitingTime = fwaitingTime        
            self.receiveFrontIndex = nResponseSliceFrontIndex
            self.receiveBackIndex = nResponseSliceBackIndex
            return 1
        except:
            return 0
 
    def Connect(self, strPort = None, nBaudrate = 9600, strParity = 'N', nDatabits = 8, nStopbits = 1, nTimeout = 1):
        if strPort is not None:
            self.com.port = strPort
        
        if nBaudrate != 9600:
            self.com.baudrate = nBaudrate
        
        if strParity != 'N':
            self.com.parity = strParity
        
        if nDatabits != 8:
            self.com.bytesize = nDatabits
        
        if nStopbits != 1:
            self.com.stopbits = nStopbits
        
        if nTimeout != 1:
            self.com.timeout = nTimeout
        
        try:
            self.com.open()
        except Exception as e:
            print("연결 실패.", e)
            return 0
        
        if self.com.isOpen():
            print("연결 성공")
            return 1
        else:
            print("연결 실패")
            return 0

    def Send(self, sendStr, suffix = None):
        if isinstance(sendStr, bytes):
            bSendStr = sendStr            
        else:
            bSendStr = sendStr.encode('ascii')

        if suffix is not None:
            self.com.write(bSendStr + suffix)
        else:
            self.com.write(bSendStr + self.sendSuffix)

        # Serial Line clear
        self.com.flushInput()
    
    def ReceiveStr(self, waitingTime = None, strLength = None, suffix = None, frontSliceIndex = None, backSliceIndex = None):
        # funcTimeStart = time.perf_counter()
        if waitingTime is None:
            waitingTime = self.waitingTime

        if strLength is None:
            strLength = self.receiveMaxLength
        
        if suffix is None:
            suffix = self.receiveSuffix
        
        if frontSliceIndex == None:
            frontSliceIndex = self.receiveFrontIndex

        if backSliceIndex == None:
            backSliceIndex = self.receiveBackIndex
        # print(time.perf_counter() - funcTimeStart)

        # 처리 시간이 긴 명령 테스트 필요
        bRecevieStr = bytearray()
        bRecevieStrLength = 0

        start = time.perf_counter()
        while bRecevieStrLength < strLength:
            bTempChar = self.com.read(1)
            bRecevieStr += bTempChar
            bRecevieStrLength += 1
            if bTempChar == suffix:
                # print("Response time : %f"%(time.perf_counter() - start))
                break
                    
            loopTime = time.perf_counter() - start
            if loopTime >= waitingTime:
                # print("waiting time out : %f"%(loopTime))
                return "timeout"

        # print(bRecevieStr)
        # return bRecevieStr[frontSliceIndex:backSliceIndex].decode('ascii')
        return bRecevieStr[frontSliceIndex:backSliceIndex].decode()

    def Disconnect(self):

        self.com.close()
        if not self.com.isOpen():
            print("종료 성공")
            return 1
        else:
            print("종료 실패")
            return 0

    def ManualCommand_T(self):

        while True:
            userInput = input("Manual Command : ")
            if userInput == 'exit':
                self.Disconnect()
                break
            else:
                self.Send(userInput)
                print(self.ReceiveStr())

    def StepCommand(self, sendStr, sendSuffix = None, waitingTime = None, strLength = None, receiveSuffix = None,
    frontSliceIndex = None, backSliceIndex = None, boolResponse = True):

        responseStr = ''
        self.Send(sendStr, sendSuffix)
        # send & receive interval time
        time.sleep(0.15)
        if boolResponse == 'True' or boolResponse == 'true' or boolResponse == True:
            responseStr = self.ReceiveStr(waitingTime, strLength, receiveSuffix, frontSliceIndex, backSliceIndex)

        return responseStr

    def LineClear(self):
        print("transmit buffer : ", self.com.out_waiting)
        print("Receive buffer : ", self.com.in_waiting)
        print("--------------Clear--------------")
        self.com.reset_input_buffer()
        self.com.reset_output_buffer()
        print("transmit buffer : ", self.com.out_waiting)
        print("Receive buffer : ", self.com.in_waiting)


class AntMono(AntSerial):
    def __init__(self, strPort=None, nBaudrate=9600, strParity='N', nDatabits=8, nStopbits=1, nTimeout=1, bXonxoff=False):
        super().__init__(strPort, nBaudrate, strParity, nDatabits, nStopbits, nTimeout, bXonxoff)

    def open(self):
        """ open """
        if self.Connect() == 1:
            self.CommunicationRule(b'\r', b'>', nResponseSliceFrontIndex=2, nResponseSliceBackIndex=-1)
            return True
        else:
            return False
    
    def close(self):
        """ close """
        if self.Disconnect() == 1:
            return True
        else:
            return False
    
    def get_version(self):
        """ get_version """
        if not self.com.is_open:
            return False
        
        self.Send("?VER")
        return self.ReceiveStr()

    def get_grating(self):
        """ get_grating """
        if not self.com.is_open:
            return False
        
        self.Send("?GRAT")
        return self.ReceiveStr()

    def set_grating(self, grating):
        """ set_grating """
        if not self.com.is_open:
            return False
        
        self.Send("!GRAT %d"%(grating))
        return self.ReceiveStr()
    
    def get_filter(self):
        """ get_filter """
        if not self.com.is_open:
            return False
        
        self.Send("?FILT1")
        return self.ReceiveStr()

    def set_filter(self, filter):
        """ set_filter """
        if not self.com.is_open:
            return False
        
        self.Send("!FILT1 %d"%(filter))
        return self.ReceiveStr()

    def get_wavelength(self):
        """ get_wavelength """
        if not self.com.is_open:
            return False
        
        self.Send("?PW")
        return self.ReceiveStr()

    def set_wavelength(self, wavelength):
        """ set_wavelength\n
        wavelength: float [nm]
        """
        if not self.com.is_open:
            return False
        
        self.Send("!GW %0.1f"%(wavelength))
        return self.ReceiveStr()

    def open_shutter(self):
        """ open_shutter """
        if not self.com.is_open:
            return False
        
        self.Send("!SHUTTER 0")
        return True

    def close_shutter(self):
        """ open_shutter """
        if not self.com.is_open:
            return False
        
        self.Send("!SHUTTER 1")
        return False

    def get_outport(self):
        """ get_outport"""
        if not self.com.is_open:
            return False
        
        self.Send("?PORTOUT")
        return self.ReceiveStr()

    def set_outport(self, port):
        """
        port: string (A, B, C)
        """
        if not self.com.is_open:
            return False
        
        self.Send("!PORTOUT %s"%(port))
        return self.ReceiveStr()
    
    def get_change_filter(self):
        """ get_change_filter """
        if not self.com.is_open:
            return False
        
        self.Send("?CHNGF1")
        return self.ReceiveStr()

    def set_change_filter(self, change_filter_query):
        """ set_change_filter """
        if not self.com.is_open:
            return False
        
        self.Send("=CHNGF1 %s"%(change_filter_query))
        return self.ReceiveStr()

    def get_change_grating(self):
        """ get_change_grating """
        if not self.com.is_open:
            return False
        
        self.Send("?CHNGGR")
        return self.ReceiveStr()

    def set_change_grating(self, change_grating_query):
        """ set_change_grating """
        if not self.com.is_open:
            return False
        
        self.Send("=CHNGGR %s"%(change_grating_query))
        return self.ReceiveStr()

    def get_bandpass(self):
        """ get_bandpass """
        if not self.com.is_open:
            return False
        
        self.Send("?BANDPASS")
        return self.ReceiveStr()

    def set_bandpass(self, bandpass):
        """ set_bandpass\n
        bandpass: float [nm]
        """
        if not self.com.is_open:
            return False
        
        self.Send("=BANDPASS %0.1f"%(bandpass))
        return self.ReceiveStr()

    def get_slit(self, port):
        """ get_slit\n
        port: string (A, B, C)
        """
        if not self.com.is_open:
            return False
        
        self.Send("?SLIT%s"%(port))
        return self.ReceiveStr()

    def set_slit(self, port, width):
        """ set_slit\n
        port: string (A, B, C)\n
        width: integer [um]
        """
        if not self.com.is_open:
            return False
        
        self.Send("!SLIT%s %d"%(port, width))
        return self.ReceiveStr()

# serial port 정보 입력
mono = AntMono('/dev/tty.usbserial-1420')
mono.open()

print(mono.get_version())
print(mono.set_outport("C"))
print(mono.get_slit("A"))
print(mono.get_slit("B"))
print(mono.get_slit("C"))

