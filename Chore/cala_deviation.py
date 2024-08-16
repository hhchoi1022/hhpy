#%%
import pandas as pd
import numpy as np
import datetime
import os
from openpyxl import load_workbook

class Export:
    
    def __init__(self,
                 idx,
                 date_bl,
                 date_inv,
                 num_inv,
                 country,
                 key_cnt,
                 price_tmp,
                 price_cert,
                 db_exchange,
                 db_export):
        self.idx = idx
        self.date_bl = date_bl
        self.date_inv = date_inv
        self.date = self.date_bl
        if pd.isnull(self.date_bl):
            self.date = self.date_inv
        self.num_inv = num_inv[-6:]
        self.country = country
        self.key_cnt = key_cnt
        self.price_tmp = price_tmp
        self.price_cert = price_cert
        self.db_exchange = db_exchange
        self.db_export = db_export
        self.price = self.price_cert
        if pd.isnull(self.price_cert):
            self.price = self.price_tmp
        self.price_remain = self.price
        self.income_list = []
        
    @property
    def exchange(self):
        date = self.date
        value = float(self.db_exchange[self.db_exchange['날짜(2년치, 매년 수정할것)'] == date]['환율'])
        return value

    @property
    def exportNum(self,
                  key_export : str ='1098'):
        try:
            val = self.db_export.loc[self.num_inv]
            value = val.values[np.flatnonzero(np.char.startswith(val.values.astype(str), key_export))][0]
        except:
            value = '국내계류중'
        return value
    
    def __repr__(self):
        return f'Export(Index:{self.idx}, Date:{str(self.date)[:10]}, Remaining:{self.price_remain}/{self.price})'
    

class Income:
    
    def __init__(self,
                 idx,
                 date_in,
                 transfer_USD,
                 transfer_charge,
                 transfer_charge_bank,
                 ID,
                 income_USD,
                 income_exchange,
                 exchange_USD,
                 exchange_exchange,
                 exchange_date,
                 income_KRW,
                 ):
        self.idx = idx
        self.date_in = date_in
        self.transfer_USD = transfer_USD
        self.transfer_charge = transfer_charge
        if pd.isnull(transfer_charge):
            self.transfer_charge = 0
        self.transfer_charge_bank = transfer_charge_bank
        if pd.isnull(transfer_charge_bank):
            self.transfer_charge_bank =0
        self.ID = ID
        self.income_USD = income_USD
        self.income_exchange = float(str(income_exchange).replace(',',''))
        self.exchange_USD = exchange_USD
        self.exchange_exchange = exchange_exchange
        self.exchange_date = exchange_date
        self.income_KRW = income_KRW
        self.price_remain = transfer_USD
        self.export_list = []
    
    def __repr__(self):
        return f'Income(Index:{self.idx}, Date:{str(self.date_in)[:10]}, Remaining:{self.price_remain}/{self.transfer_USD})'
    
class Section(object):
    pass

class Database:
    
    def __init__(self,
                 file_name = './출고내역정리_22.xlsx',
                 **kwargs):
        self.file_name = file_name
        self.db = self._load_DB()
        self.all_sheet_names = self._get_country_list()
        self.db_export = self._load_DB_export()
        self.db_exchange = self._load_DB_exchange()

    def _load_DB(self):
        db = pd.read_excel(self.file_name, engine = 'openpyxl', sheet_name=None)
        return db
    
    def _load_DB_export(self):
        data = self.db['수출신고']
        keys = data.T.iloc[0]
        db_export = data.copy()
        db_export.index = keys
        return db_export
    
    def _load_DB_exchange(self):
        data = self.db['환율']
        keys = data['날짜(2년치, 매년 수정할것)']
        db_exchange = data.copy()
        db_exchange.index = keys
        return db_exchange
    
    def _get_country_list(self):
        return list(['CKH,LKB', 'NGM', 'INA', 'INT,IRR', 'PMR', 'PKK', 'RK,RKH', 'UMH,UQS', 'ZBJ'])
    '''    
    def update(self,
               savepath : str = './history/'):
        if not os.path.isdir(savepath):
            os.makedirs(savepath)
        all_db = self.db
        for country in self.all_sheet_names:
            sheet = Sheet(sheet_name = country)
            all_db[country] = sheet.data.data_updated
        history_filename = f"{savepath}{datetime.datetime.now().strftime('%y%m%d_%H%M%S')}_"+os.path.basename(self.file_name)
        os.system(f'cp {self.file_name} {history_filename}')
        writer = pd.ExcelWriter(self.file_name, engine = 'openpyxl')
        for key, values in all_db.items():
            values.to_excel(writer, sheet_name = key, index = False)
        writer.save()
        return all_db
    '''
    def update(self,
               savepath : str = './history/'):
        if not os.path.isdir(savepath):
            os.makedirs(savepath)
        
        all_db = self.db
        for country in self.all_sheet_names:
            sheet = Sheet(sheet_name = country)
            all_db[country] = sheet.data.data_updated
        
        book = load_workbook(self.file_name)
        writer = pd.ExcelWriter(self.file_name, engine='openpyxl')
        writer.book = book
        history_filename = f"{savepath}{datetime.datetime.now().strftime('%y%m%d_%H%M%S')}_"+os.path.basename(self.file_name)
        os.system(f'cp {self.file_name} {history_filename}')
        for sheet in writer.book.worksheets:
            writer.book.remove(sheet)
        # Write the updated DataFrames back to the Excel file
        for sheet_name, df in all_db.items():
            df.to_excel(writer, index=False, sheet_name=sheet_name)
        writer.save()

        
#%%        
        
        
class Sheet(Database):
    
    def __init__(self,
                 sheet_name = 'NGM',
                 **kwargs):
        super().__init__()
        self.sheet_name = sheet_name
        self.data = Section()
        if not sheet_name is None:
            self.data.data, self.data.data_raw, self.data._row_1st = self.load_sheet()
            self.data.data_result = self.get_result_data()
            self.data.data_export = self.get_export_data()
            self.data.data_income = self.get_income_data()
            self.data.data_export_detail = self.get_export_data_detail()
            self.data.data_income_detail = self.get_income_data_detail()
            self.match_export_income()
            self.data.data_updated = self.update_data()
            
            
    def load_sheet(self):
        data = self.db[self.sheet_name]
        row_1st = data.columns
        header = data.iloc[0]
        df = data[1:]
        df.rename(columns = header, inplace = True)
        return df, data, row_1st
    
    def get_export_data_detail(self,
                              idx_min : int = 0,
                              idx_max : int = 7):
        db_sheet = self.data.data
        data = db_sheet.T.iloc[idx_min:idx_max].T
        idx = data['INV.dat'].dropna().index
        export_data = data.loc[idx]
        all_export = []
        for idx, row in export_data.iterrows():
            arg = [idx] + row.values.tolist()
            income = Export(*arg, db_exchange= self.db_exchange, db_export= self.db_export)
            all_export.append(income)
        return all_export
    
    def get_income_data_detail(self,
                              idx_min : int = 9,
                              idx_max : int = 20):
        db_sheet = self.data.data
        data = db_sheet.T.iloc[idx_min:idx_max].T
        idx = data['입금일'].dropna().index
        income_data = data.loc[idx]
        all_income = []
        for idx, row in income_data.iterrows():
            arg = [idx] + row.values.tolist()
            income = Income(*arg)
            all_income.append(income)
        return all_income

    def get_export_data(self,
                        idx_min : int = 0,
                        idx_max : int = 7):
        db_sheet = self.data.data
        data = db_sheet.T.iloc[idx_min:idx_max].T
        idx = data['INV.dat'].dropna().index
        return data.loc[idx]
        
    def get_income_data(self,
                        idx_min : int = 9,
                        idx_max : int = 20):
        db_sheet = self.data.data
        data = db_sheet.T.iloc[idx_min:idx_max].T
        idx = data['입금일'].dropna().index
        return data.loc[idx]

    def get_result_data(self,
                        idx_min : int = 22,
                        idx_max : int = 33 ):
        db_sheet = self.data.data
        data = db_sheet.T.iloc[idx_min:idx_max]
        result_data = data.dropna(axis = 'columns').T
        return result_data
    
    def update_data(self,
                    idx_min : int = 22,
                    idx_max : int = 33):
        db_sheet = self.data.data
        data = db_sheet.T.iloc[idx_min:idx_max].T
        update_table = self.update_result_data()
        for i, row in update_table.iterrows():
            data.loc[i] = row
        db_sheet_transpose = db_sheet.T
        for i,j in enumerate(np.arange(idx_min, idx_max)):
            db_sheet_transpose.iloc[j] = data.T.iloc[i]
        updated_db = db_sheet_transpose.T
        updated_db.loc[0] = list(db_sheet.columns)
        sorted_db = updated_db.sort_index()
        sorted_db.columns = self.data.data_raw.columns
        return sorted_db
    
    def fill_price_USD(self):
        USDpricelist = []
        for price_tmp, price_cert in zip(self.data.data_ship['임시단가'].values, self.data.data_ship['필증단가'].values):
            price = price_cert
            if np.isnan(price_cert):
                price = price_tmp
            USDpricelist.append(price)
        return USDpricelist
    
    def match_export_income(self):
        for export in self.data.data_export_detail:
            all_incomes = [income for income in self.data.data_income_detail if income.price_remain > 0]
            for income in all_incomes:
                if export.price_remain == 0:
                    break  
                price_difference = export.price_remain - income.price_remain
                if price_difference > 0:
                    export.income_list.append([income, export.price_remain-price_difference])
                    income.export_list.append([export, export.price_remain-price_difference])
                    export.price_remain -= income.price_remain
                    income.price_remain -= income.price_remain
                elif price_difference <0:
                    export.income_list.append([income, income.price_remain+price_difference])
                    income.export_list.append([export, income.price_remain+price_difference])
                    income.price_remain -= export.price_remain
                    export.price_remain -= export.price_remain
                else:
                    export.income_list.append([income, export.price_remain-price_difference])
                    income.export_list.append([export, export.price_remain-price_difference])
                    export.price_remain -= income.price_remain
                    income.price_remain -= income.price_remain
        
    def update_result_data(self):
        columns = self.data.data_result.columns
        result_tbl = pd.DataFrame(columns = columns)
        result_tbl['BL.dat'] = self.data.data_export['BL.dat']
        for export in self.data.data_export_detail:
            prepaid = np.sum([income[0].price_remain for income in export.income_list])
            index = export.idx
            result_tbl.loc[index,'송장번호'] = export.num_inv
            result_tbl.loc[index,'신고번호'] = export.exportNum
            result_tbl.loc[index,'환율(B/L)'] = export.exchange
            result_tbl.loc[index,'단가(USD)'] = export.price
            if prepaid >0:
                result_tbl.loc[index,'입금(USD)'] = export.price - export.price_remain + prepaid
            else:
                result_tbl.loc[index,'입금(USD)'] = export.price - export.price_remain
            result_tbl.loc[index,'미수(USD)'] = export.price_remain - prepaid
            result_tbl.loc[index,'단가(KRW)'] = export.exchange * export.price
            summation = 0
            if not len(export.income_list) == 0:
                for income, amount_income in export.income_list:
                    price = (amount_income - income.transfer_charge - income.transfer_charge_bank) + income.price_remain
                    summation += price*income.income_exchange
            result_tbl.loc[index,'입금(KRW)'] = summation
            result_tbl.loc[index,'미수(KRW)'] = export.exchange * (export.price_remain - prepaid)
            result_tbl.loc[index,'환차손익'] = summation - export.exchange * export.price + export.exchange * (export.price_remain - prepaid)
        return result_tbl
#%%
db = Database()
A = db.update()


# %%
