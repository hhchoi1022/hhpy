class Folder:
    def __init__(self, 
                 watchdir = '/data6/obsdata/KCT_STX16803/zip', 
                 cpdir = '/data6/obsdata/KCT_STX16803',
                 check_period = 10,
                 recipients = ['hhchoi1022@gmail.com']):
        
        self.watchdir = watchdir
        self.cpdir = cpdir
        self.tmpdir = f'{watchdir}/temp'
        os.makedirs(self.tmpdir, exist_ok= True)
    
    def change_hdr(self, folderpath):
        filelist = glob.glob(f"{folderpath}/*")
        for filepath in filelist:            
            try:
                hdr = fits.getheader(filepath)
                if 'OBJECT' not in hdr.keys():
                    fits.setval(filepath, 'OBJECT', value = '')
                    print(f'{os.path.basename(filepath)} is modified')
            except:
                pass
    
    def unzip(self, event_path, modify_header : bool = True):
        import datetime
        os.system(f'chmod 777 {event_path}')
        os.system(f'unzip -d {self.tmpdir} {event_path}')
        foldername = os.path.basename(event_path).split('.')[0]
        folderpath = f'{self.tmpdir}/{foldername}'
        if modify_header:
            self.change_hdr(folderpath)
        os.system(f'mv {folderpath} {self.cpdir}')
        os.system(f'chmod 777 {folderpath}')
        