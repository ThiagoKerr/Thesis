# -*- coding: utf-8 -*-
"""
Processador Geodésico - Gravidade + Grade Geoidal (Interpolação Bilinear)

Entrada 1 (gravimetria): TXT/CSV/XLSX com colunas (mapeáveis via GUI):
    - ID
    - lat_tide_free (graus)
    - lon_tide_free (graus)
    - h_tide_free (m)
    - g_mean_tide (mGal)
    - fonte

Entrada 2 (modelo geoidal): TXT/CSV com três colunas: longitude, latitude, N (m).
    - Grade regular (p.ex., passo de 5'), cobrindo a área de interesse.

Saída XLSX com:
    - N interpolado (bilinear) a partir da grade
    - H_zero_tide = h_mean_tide - N, H_mean_tide (= H_zero_tide por padrão)
    - Demais grandezas: γ0, γ(H), correção atmosférica, etc.
"""

import os
import math
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Tuple, Dict
import numpy as np
import pandas as pd


# =============================
# CONSTANTES
# =============================
A_WGS84 = 6378137.0
F_INV   = 298.257223563
F       = 1.0 / F_INV
B_WGS84 = A_WGS84 * (1.0 - F)
E2      = F * (2.0 - F)

G_E = 9.7803253359
G_P = 9.8321849378
K   = (B_WGS84 * G_P - A_WGS84 * G_E) / (A_WGS84 * G_E)

MS2_TO_MGAL = 1e5

ATM_DELTA0_MGAL = 0.87
ATM_SCALE_M     = 8435.0

FA_GRAD_1 = 0.3087691
FA_GRAD_2 = 0.000000439

PT_C0_MGAL = -0.3086

DEFAULT_APPLY_PT_ON_G    = True
DEFAULT_APPLY_PT_ON_GAM  = True
DEFAULT_H_MT_EQUALS_H_TF = True
DEFAULT_OOB_POLICY       = 'clamp'

# =============================
# FUNÇÕES
# =============================
def somigliana_gamma(phi_rad: float) -> float:
    s = math.sin(phi_rad)
    denom = math.sqrt(1.0 - E2 * s * s)
    return G_E * (1.0 + K * s * s) / denom

def geodetic_to_ecef(lat_deg, lon_deg, h_m):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    s = math.sin(lat); c = math.cos(lat)
    N = A_WGS84 / math.sqrt(1.0 - E2 * s * s)
    X = (N + h_m) * c * math.cos(lon)
    Y = (N + h_m) * c * math.sin(lon)
    Z = (N * (1.0 - E2) + h_m) * s
    return X, Y, Z

def geocentric_latitude(lat_deg, lon_deg, h_m):
    X, Y, Z = geodetic_to_ecef(lat_deg, lon_deg, h_m)
    return math.degrees(math.atan2(Z, math.hypot(X, Y)))

def atmospheric_correction_mgal(H_m: float) -> float:
    return ATM_DELTA0_MGAL * math.exp(-max(0.0, H_m) / ATM_SCALE_M)

def g_pt_delta_mgal(phi_rad: float) -> float:
    return PT_C0_MGAL * (1.0 - 1.5 * math.sin(phi_rad)**2)

def normal_gravity_on_surface_mgal(phi_rad: float, H_m: float) -> float:
    gamma0_ms2 = somigliana_gamma(phi_rad)
    gamma0_mgal = gamma0_ms2 * MS2_TO_MGAL
    return gamma0_mgal - FA_GRAD_1 * H_m + FA_GRAD_2 * (H_m**2)

def read_table_auto_with_header(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, engine='python')
    except Exception:
        try:
            return pd.read_excel(path)
        except Exception:
            df = pd.read_csv(path, engine='python', header=None)
            df.columns = [f'col{i}' for i in range(df.shape[1])]
            return df

# =============================
# GEOID GRID
# =============================
def read_geoid_grid_lon_lat_N(path: str):
    """
    Lê SAM_GEOID (ou equivalente) com 3 colunas: lon lat N (graus decimais, N em m).
    Não altera nem arredonda valores. Apenas organiza em uma grade (nlat x nlon).
    """
    import pandas as pd
    import numpy as np

    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python", names=["lon","lat","N"])

    # Conjuntos únicos ordenados (sem arredondamentos)
    lons = np.sort(df["lon"].unique())
    lats = np.sort(df["lat"].unique())

    # Monta matriz N(lat,lon)
    nlat, nlon = len(lats), len(lons)
    grid = np.full((nlat, nlon), np.nan, dtype=float)
    lon_to_ix = {v: i for i, v in enumerate(lons)}
    lat_to_ix = {v: i for i, v in enumerate(lats)}

    for lon, lat, val in df.itertuples(index=False):
        i = lat_to_ix[lat]
        j = lon_to_ix[lon]
        grid[i, j] = float(val)

    # Se houver NaN, **não preenche**: isso indica falta de par (lon,lat) no arquivo.
    if np.isnan(grid).any():
        raise ValueError("A grade geoidal possui lacunas (pares lon/lat ausentes). Verifique o arquivo.")

    return lons, lats, grid



def bilinear_interpolate_grid(lons, lats, grid, qlon, qlat, oob_policy='nan'):
    import numpy as np

    qlon = np.asarray(qlon, dtype=float)
    qlat = np.asarray(qlat, dtype=float)

    ix = np.searchsorted(lons, qlon) - 1
    iy = np.searchsorted(lats, qlat) - 1

    mask_oob = (qlon < lons[0]) | (qlon > lons[-1]) | (qlat < lats[0]) | (qlat > lats[-1])

    ix = np.clip(ix, 0, len(lons) - 2)
    iy = np.clip(iy, 0, len(lats) - 2)

    x0, x1 = lons[ix],     lons[ix + 1]
    y0, y1 = lats[iy],     lats[iy + 1]

    tx = np.where(x1 != x0, (qlon - x0) / (x1 - x0), 0.0)
    ty = np.where(y1 != y0, (qlat - y0) / (y1 - y0), 0.0)

    f00 = grid[iy,     ix]
    f10 = grid[iy,     ix + 1]
    f01 = grid[iy + 1, ix]
    f11 = grid[iy + 1, ix + 1]

    a = f00 * (1 - tx) + f10 * tx
    b = f01 * (1 - tx) + f11 * tx
    out = a * (1 - ty) + b * ty

    if oob_policy == 'nan':
        out = out.astype(float)
        out[mask_oob] = np.nan
    return out


def normalize_lon_to_grid(lon_array, grid_lons):
    """
    Ajusta somente as longitudes de entrada (PPTE) para o mesmo domínio das longitudes da grade,
    sem alterar a grade. Se a grade está em 0–360, projeta PPTE para 0–360; se está em −180..180,
    projeta PPTE para −180..180.
    """
    import numpy as np
    lon = np.asarray(lon_array, dtype=float)
    glon_min, glon_max = float(grid_lons.min()), float(grid_lons.max())
    # Grade em 0–360?
    if glon_min >= 0.0 and glon_max <= 360.0:
        lon = lon % 360.0
        lon = np.where(lon < 0.0, lon + 360.0, lon)
    else:
        lon = ((lon + 180.0) % 360.0) - 180.0
    return lon

# =============================
# PROCESSAMENTO
# =============================
def process_with_grid(grav_df, lons, lats, gridN, opts):
    df = grav_df.copy()
    # alinhar domínio de longitude do PPTE ao domínio da grade (sem tocar na grade!)
    df = grav_df.copy()
    df['lon_tide_free'] = normalize_lon_to_grid(df['lon_tide_free'].values, lons)

    Nvals = bilinear_interpolate_grid(
        lons, lats, gridN,
        df['lon_tide_free'].values.astype(float),
        df['lat_tide_free'].values.astype(float),
        oob_policy='nan'
    )
    df['N'] = Nvals
    if np.isnan(df['N']).any():
        raise ValueError(
            "Há pontos fora da extensão da grade geoidal.\n"
            f"Lat PPTE: [{df['lat_tide_free'].min():.6f}, {df['lat_tide_free'].max():.6f}]  "
            f"vs Grade: [{lats.min():.6f}, {lats.max():.6f}]\n"
            f"Lon PPTE (ajustada ao domínio): [{df['lon_tide_free'].min():.6f}, {df['lon_tide_free'].max():.6f}]  "
            f"vs Grade: [{lons.min():.6f}, {lons.max():.6f}]"
        )

    # Interpola N
    Nvals = bilinear_interpolate_grid(
        lons, lats, gridN,
        df['lon_tide_free'].values.astype(float),
        df['lat_tide_free'].values.astype(float),
        oob_policy='nan'  # <- fundamental!
    )
    df['N'] = Nvals

    # Se houver NaN, aponta extensão e aborta (evita resultado silenciosamente errado)
    if np.isnan(df['N']).any():
        raise ValueError(
            "Há pontos fora da extensão da grade geoidal.\n"
            f"Lat PPTE: [{df['lat_tide_free'].min():.6f}, {df['lat_tide_free'].max():.6f}]  "
            f"vs Grade: [{lats.min():.6f}, {lats.max():.6f}]\n"
            f"Lon PPTE: [{df['lon_tide_free'].min():.6f}, {df['lon_tide_free'].max():.6f}]  "
            f"vs Grade: [{lons.min():.6f}, {lons.max():.6f}]\n"
            "Carregue uma grade que cubra todos os pontos."
        )

    # Demais cálculos (como você já tinha)
    df['lat_geocentrica_deg'] = [
        geocentric_latitude(lat, lon, h) for lat,lon,h in
        zip(df['lat_tide_free'], df['lon_tide_free'], df['h_tide_free'])
    ]
    df['h_mean_tide'] = df['h_tide_free']
    df['H_zero_tide'] = df['h_mean_tide'] - df['N']
    df['H_mean_tide'] = df['H_zero_tide']
    df['delta_g_atm_mt_mGal'] = df['H_mean_tide'].apply(atmospheric_correction_mgal)
    df['g_atm_mean_tide_mGal'] = df['g_mean_tide'] + df['delta_g_atm_mt_mGal']
    phi_rad = np.radians(df['lat_tide_free'].astype(float))
    delta_pt = np.array([g_pt_delta_mgal(p) for p in phi_rad])
    df['g_atm_zero_tide_mGal'] = df['g_atm_mean_tide_mGal'] - delta_pt
    df['gamma0_mGal'] = [somigliana_gamma(math.radians(phi))*MS2_TO_MGAL for phi in df['lat_tide_free']]
    df['gamma_surface_mean_tide_mGal'] = [
        normal_gravity_on_surface_mgal(math.radians(phi), H)
        for phi,H in zip(df['lat_tide_free'], df['H_mean_tide'])
    ]
    df['gamma_surface_zero_tide_mGal'] = df['gamma_surface_mean_tide_mGal'] - delta_pt
    df['disturbio_g_zero_tide_mGal'] = df['g_atm_zero_tide_mGal'] - df['gamma_surface_zero_tide_mGal']
    rename_map = {
        'lat_tide_free': 'Geodetic Lat (Tide Free)',
        'lon_tide_free': 'Geodetic Lon (Tide Free)',
        'h_tide_free': 'h (Tide Free)',
        'g_mean_tide': 'g (mean tide)',
        'N': 'N (Zero Tide)',
        'lat_geocentrica_deg': 'Geocentric Lat (Tide Free)',
        'h_mean_tide': 'h (Mean Tide)',
        'H_zero_tide': 'H (Zero Tide)',
        'H_mean_tide': 'H (Mean Tide)',
        'delta_g_atm_mt_mGal': 'Atm Correction (Mean Tide)',
        'g_atm_mean_tide_mGal': 'g with atm correction (Mean Tide)',
        'g_atm_zero_tide_mGal': 'g with atm correction (Zero Tide)',
        'gamma0_mGal': 'Normal gravity on the ellipsoid (Tide Free).',
        'gamma_surface_mean_tide_mGal': 'Normal Gravity on the surface (Mean Tide)',
        'gamma_surface_zero_tide_mGal': 'Normal gravity on the surface (Zero Tide).',
        'disturbio_g_zero_tide_mGal': 'Gravity disturbance (Zero Tide).'
    }
    df = df.rename(columns=rename_map)
    return df



# =============================
# GUI
# =============================
def read_grav_file(path: str) -> pd.DataFrame:
    # Lê Excel/CSV; usa as 5 primeiras colunas nesta ordem: ID, lat, lon, h, g
    try:
        df = pd.read_excel(path, header=0)
    except Exception:
        df = pd.read_csv(path, engine='python')

    df = df.iloc[:, :5].copy()
    df.columns = ['ID', 'lat_tide_free', 'lon_tide_free', 'h_tide_free', 'g_mean_tide']
    df['fonte'] = ""

    # normaliza tipos
    for c in ['lat_tide_free', 'lon_tide_free', 'h_tide_free', 'g_mean_tide']:
        df[c] = pd.to_numeric(df[c], errors='coerce')

    # normaliza longitude para −180..180 e garante Oeste negativo
    # (se vier 0–360, converte; se vier +W por engano, troca sinal)
    lon = df['lon_tide_free'].values.astype(float)
    lon = np.where(lon > 180.0, lon - 360.0, lon)   # 0–360 -> −180..180
    # Heurística: se >90% são positivos em região que deveria ser W, inverte sinal
    if (lon > 0).mean() > 0.9:
        lon = -np.abs(lon)
    df['lon_tide_free'] = lon

    return df

def write_named_parameters(workbook, worksheet):
    # escreve cabeçalho e valores na aba 'parametros' e cria nomes definidos
    params = [
        ("A_WGS84", A_WGS84),
        ("F", F),
        ("E2", E2),
        ("G_E", G_E),
        ("G_P", G_P),
        ("K", K),
        ("MS2_TO_MGAL", MS2_TO_MGAL),
        ("FA_GRAD_1", FA_GRAD_1),
        ("FA_GRAD_2", FA_GRAD_2),
        ("ATM_DELTA0_MGAL", ATM_DELTA0_MGAL),
        ("ATM_SCALE_M", ATM_SCALE_M),
        ("PT_C0_MGAL", PT_C0_MGAL),
    ]
    worksheet.write(0, 0, "Parametro")
    worksheet.write(0, 1, "Valor")
    for i, (name, val) in enumerate(params, start=1):
        worksheet.write(i, 0, name)
        worksheet.write_number(i, 1, float(val))
        # define nome global apontando para a célula de valor (linha i+1 em Excel 1-based)
        workbook.define_name(f"{name}", f"=parametros!$B${i+1}")


def write_with_constants_sheet(out_df, writer):
    """
    Cria:
      - 'resultado' com dados + fórmulas
      - 'constantes' com parâmetros numéricos visíveis
    """
    wb = writer.book
    ws_res = wb.add_worksheet("resultado")
    ws_cst = wb.add_worksheet("constantes")

    # ==========================
    # 1. escreve constantes na aba 'constantes'
    # ==========================
    params = [
        ("A_WGS84", A_WGS84),
        ("E2", E2),
        ("G_E", G_E),
        ("K", K),
        ("MS2_TO_MGAL", MS2_TO_MGAL),
        ("FA_GRAD_1", FA_GRAD_1),
        ("FA_GRAD_2", FA_GRAD_2),
        ("ATM_DELTA0_MGAL", ATM_DELTA0_MGAL),
        ("ATM_SCALE_M", ATM_SCALE_M),
        ("PT_C0_MGAL", PT_C0_MGAL),
    ]
    ws_cst.write(0, 0, "Parametro")
    ws_cst.write(0, 1, "Valor")
    for i, (name, val) in enumerate(params, start=1):
        ws_cst.write(i, 0, name)
        ws_cst.write_number(i, 1, float(val))

    # Função auxiliar: retorna referência absoluta na aba constantes
    def ref_const(row):
        return f"constantes!$B${row}"

    # refs
    A_WGS84_c   = ref_const(1+1)   # linha 2
    E2_c        = ref_const(2+1)   # linha 3
    G_E_c       = ref_const(3+1)
    K_c         = ref_const(4+1)
    MS2_c       = ref_const(5+1)
    FA1_c       = ref_const(6+1)
    FA2_c       = ref_const(7+1)
    ATM0_c      = ref_const(8+1)
    ATMscale_c  = ref_const(9+1)
    PT_c        = ref_const(10+1)

    # ==========================
    # 2. cabeçalhos da aba resultado
    # ==========================
    headers = [
        "ID","Geodetic Lat (Tide Free)","Geodetic Lon (Tide Free)",
        "h (Tide Free)","g (mean tide)","N (Zero Tide)",
        "Geocentric Lat (Tide Free)",
        "h (Mean Tide)","H (Zero Tide)","H (Mean Tide)",
        "Atm Correction (Mean Tide)","g with atm correction (Mean Tide)",
        "delta_g_PT (mGal)","g with atm correction (Zero Tide)",
        "Normal gravity on the ellipsoid (Tide Free).",
        "Normal Gravity on the surface (Mean Tide)",
        "Normal gravity on the surface (Zero Tide).",
        "Gravity disturbance (Zero Tide)."
    ]
    for j, h in enumerate(headers):
        ws_res.write(0, j, h)

    # ==========================
    # 3. escreve dados + fórmulas
    # ==========================
    cols_in = [
        "ID",
        "Geodetic Lat (Tide Free)",
        "Geodetic Lon (Tide Free)",
        "h (Tide Free)",
        "g (mean tide)",
        "N (Zero Tide)",
    ]
    missing = [c for c in cols_in if c not in out_df.columns]
    if missing:
        raise KeyError(f"Faltam colunas na saída: {missing}")

    n = len(out_df)
    for i in range(n):
        r = i + 1  # linha em Excel
        ws_res.write(r, 0, out_df.iloc[i][cols_in[0]])
        ws_res.write_number(r, 1, float(out_df.iloc[i][cols_in[1]]))
        ws_res.write_number(r, 2, float(out_df.iloc[i][cols_in[2]]))
        ws_res.write_number(r, 3, float(out_df.iloc[i][cols_in[3]]))
        ws_res.write_number(r, 4, float(out_df.iloc[i][cols_in[4]]))
        ws_res.write_number(r, 5, float(out_df.iloc[i][cols_in[5]]))

    # Fórmulas
    for i in range(2, n+2):
        lat = f"$B${i}"
        lon = f"$C${i}"
        h   = f"$D${i}"
        g   = f"$E${i}"
        N   = f"$F${i}"

        Nphi = f"{A_WGS84_c}/SQRT(1-{E2_c}*SIN(RADIANS({lat}))^2)"
        X = f"({Nphi}+{h})*COS(RADIANS({lat}))*COS(RADIANS({lon}))"
        Y = f"({Nphi}+{h})*COS(RADIANS({lat}))*SIN(RADIANS({lon}))"
        Z = f"({Nphi}*(1-{E2_c})+{h})*SIN(RADIANS({lat}))"

        ws_res.write_formula(i-1, 6,  f"=DEGREES(ATAN2({Z},SQRT(({X})^2+({Y})^2)))")
        ws_res.write_formula(i-1, 7,  f"={h}")               # h_mean
        ws_res.write_formula(i-1, 8,  f"=D{i}-F{i}")         # H_zero
        ws_res.write_formula(i-1, 9,  f"=I{i}")              # H_mean
        ws_res.write_formula(i-1,10, f"={ATM0_c}*EXP(-MAX(0,J{i})/{ATMscale_c})")
        ws_res.write_formula(i-1,11, f"=E{i}+K{i}")
        ws_res.write_formula(i-1,12, f"={PT_c}*(1-1.5*SIN(RADIANS({lat}))^2)")
        ws_res.write_formula(i-1,13, f"=L{i}-M{i}")
        ws_res.write_formula(i-1,14, f"={G_E_c}*(1+{K_c}*SIN(RADIANS({lat}))^2)/"
                                     f"SQRT(1-{E2_c}*SIN(RADIANS({lat}))^2)*{MS2_c}")
        ws_res.write_formula(i-1,15, f"=O{i}-{FA1_c}*J{i}+{FA2_c}*J{i}^2")
        ws_res.write_formula(i-1,16, f"=P{i}-M{i}")
        ws_res.write_formula(i-1,17, f"=N{i}-Q{i}")

    ws_res.freeze_panes(1, 1)
    ws_res.autofilter(0, 0, n, len(headers)-1)





class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Processador Gravimetria + Geóide")
        self.geometry("720x500")
        self.grav_path = tk.StringVar()
        self.grid_path = tk.StringVar()
        self.out_path  = tk.StringVar(value="GravityDisturbance_Output.xlsx")
        self.colmap = {
            'ID': tk.StringVar(value='ID'),
            'lat_tide_free': tk.StringVar(value='lat_tide_free'),
            'lon_tide_free': tk.StringVar(value='lon_tide_free'),
            'h_tide_free': tk.StringVar(value='h_tide_free'),
            'g_mean_tide': tk.StringVar(value='g_mean_tide'),
            'fonte': tk.StringVar(value='fonte'),
        }
        self._build_ui()

    def _build_ui(self):
        pad={'padx':5,'pady':5}
        frm=ttk.Frame(self); frm.pack(fill='both',expand=True,**pad)
        ttk.Label(frm,text="Arquivo Gravimetria").grid(row=0,column=0,sticky='w')
        ttk.Entry(frm,textvariable=self.grav_path,width=60).grid(row=1,column=0,sticky='we')
        ttk.Button(frm,text="Procurar",command=self.browse_grav).grid(row=1,column=1)
        ttk.Label(frm,text="Grade Geoidal (lon lat N)").grid(row=2,column=0,sticky='w')
        ttk.Entry(frm,textvariable=self.grid_path,width=60).grid(row=3,column=0,sticky='we')
        ttk.Button(frm,text="Procurar",command=self.browse_grid).grid(row=3,column=1)
        ttk.Label(frm,text="Saída XLSX").grid(row=4,column=0,sticky='w')
        ttk.Entry(frm,textvariable=self.out_path,width=60).grid(row=5,column=0,sticky='we')
        ttk.Button(frm,text="Salvar como",command=self.browse_out).grid(row=5,column=1)
        ttk.Button(frm,text="Gerar XLSX",command=self.run).grid(row=6,column=1,sticky='e')

    def browse_grav(self):
        p=filedialog.askopenfilename(); 
        if p: self.grav_path.set(p)
    def browse_grid(self):
        p=filedialog.askopenfilename(); 
        if p: self.grid_path.set(p)
    def browse_out(self):
        p=filedialog.asksaveasfilename(defaultextension=".xlsx"); 
        if p: self.out_path.set(p)

    def run(self):
        try:
            grav_df = read_grav_file(self.grav_path.get())
            lons, lats, gridN = read_geoid_grid_lon_lat_N(self.grid_path.get())

            # Validação rápida
            lat_ok = (grav_df['lat_tide_free'].min() >= lats.min()-1e-9) and (grav_df['lat_tide_free'].max() <= lats.max()+1e-9)
            lon_ok = (grav_df['lon_tide_free'].min() >= lons.min()-1e-9) and (grav_df['lon_tide_free'].max() <= lons.max()+1e-9)
            if not (lat_ok and lon_ok):
                raise ValueError(
                    "Extensão da grade não cobre todos os pontos.\n"
                    f"Lat PPTE: [{grav_df['lat_tide_free'].min():.6f}, {grav_df['lat_tide_free'].max():.6f}]  "
                    f"vs Grade: [{lats.min():.6f}, {lats.max():.6f}]\n"
                    f"Lon PPTE: [{grav_df['lon_tide_free'].min():.6f}, {grav_df['lon_tide_free'].max():.6f}]  "
                    f"vs Grade: [{lons.min():.6f}, {lons.max():.6f}]"
                )

            out_df = process_with_grid(grav_df, lons, lats, gridN, opts={})
            with pd.ExcelWriter(self.out_path.get(), engine='xlsxwriter') as writer:
                write_with_constants_sheet(out_df, writer)





            messagebox.showinfo("OK", f"Arquivo gerado:\n{self.out_path.get()}")
        except Exception as e:
            messagebox.showerror("Erro", str(e))

def main():
    app=App(); app.mainloop()

if __name__=="__main__":
    main()
