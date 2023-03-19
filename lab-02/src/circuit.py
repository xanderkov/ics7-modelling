class Circuit:
    
    def __init__(self):
        self.R = 0.35
        self.I_e = 12
        self.L_k = 187 * 10**(-6)
        self.C_k = 268 * 10**(-6)
        self.R_k = 0.25
        self.U_co = 1400
        self.I_o = 0.3
        self.T_w = 2000
        
        self.table1_filename = "./data/table1.csv"
        self.table2_filename = "./data/table2.csv"