# -*- coding: utf-8 -*-
"""
Processador Geodésico - Gravidade + Grade Geoidal (Interpolação Bilinear)

Entrada 1 (gravimetria): TXT/CSV/XLSX com colunas:
    - ID
    - lat_tide_free (graus)
    - lon_tide_free (graus)
    - h_tide_free (m)
    - g_mean_tide (mGal)
    - fonte

Entrada 2 (modelo geoidal): TXT/CSV com três colunas: longitude, latitude, N (m).
    - Grade regular (p.ex., passo de 5'), cobrindo a área de interesse.

Saída XLSX:
    - N interpolado (bilinear) a partir da grade
    - Restante das colunas preenchidas com fórmulas Excel (literais)
"""

import math
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import numpy as np
import pandas as pd


# =============================
# GEOID GRID
# =============================
def read_geoid_grid_lon_lat_N(path: str):
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python", names=["lon", "lat", "N"])

    lons = np.sort(df["lon"].unique())
    lats = np.sort(df["lat"].unique())

    nlat, nlon = len(lats), len(lons)
    grid = np.full((nlat, nlon), np.nan, dtype=float)
    lon_to_ix = {v: i for i, v in enumerate(lons)}
    lat_to_ix = {v: i for i, v in enumerate(lats)}

    for lon, lat, val in df.itertuples(index=False):
        i = lat_to_ix[lat]
        j = lon_to_ix[lon]
        grid[i, j] = float(val)

    if np.isnan(grid).any():
        raise ValueError("A grade geoidal possui lacunas (pares lon/lat ausentes).")

    return lons, lats, grid


def bilinear_interpolate_grid(lons, lats, grid, qlon, qlat, oob_policy='nan'):
    qlon = np.asarray(qlon, dtype=float)
    qlat = np.asarray(qlat, dtype=float)

    ix = np.searchsorted(lons, qlon) - 1
    iy = np.searchsorted(lats, qlat) - 1

    mask_oob = (qlon < lons[0]) | (qlon > lons[-1]) | (qlat < lats[0]) | (qlat > lats[-1])

    ix = np.clip(ix, 0, len(lons) - 2)
    iy = np.clip(iy, 0, len(lats) - 2)

    x0, x1 = lons[ix], lons[ix + 1]
    y0, y1 = lats[iy], lats[iy + 1]

    tx = np.where(x1 != x0, (qlon - x0) / (x1 - x0), 0.0)
    ty = np.where(y1 != y0, (qlat - y0) / (y1 - y0), 0.0)

    f00 = grid[iy, ix]
    f10 = grid[iy, ix + 1]
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
    lon = np.asarray(lon_array, dtype=float)
    glon_min, glon_max = float(grid_lons.min()), float(grid_lons.max())
    if glon_min >= 0.0 and glon_max <= 360.0:
        lon = lon % 360.0
        lon = np.where(lon < 0.0, lon + 360.0, lon)
    else:
        lon = ((lon + 180.0) % 360.0) - 180.0
    return lon


# =============================
# PROCESSAMENTO
# =============================
def process_with_grid(grav_df, lons, lats, gridN):
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
        raise ValueError("Há pontos fora da extensão da grade geoidal.")

    rename_map = {
        'lat_tide_free': 'Geodetic Lat (Tide Free)',
        'lon_tide_free': 'Geodetic Lon (Tide Free)',
        'h_tide_free': 'h (Tide Free)',
        'g_mean_tide': 'g (mean tide)',
        'N': 'N (Zero Tide)',
    }
    df = df.rename(columns=rename_map)
    return df


# =============================
# EXPORTAÇÃO PARA EXCEL
# =============================
def write_with_constants_sheet(out_df, writer):
    """
    Cria:
      - 'resultado' com dados + FÓRMULAS (literais, em inglês: SIN/PI/ATAN/TAN).
      - 'constantes' apenas para visualização (não é referenciada nas fórmulas).
    """
    wb = writer.book
    ws_res = wb.add_worksheet("resultado")
    ws_cst = wb.add_worksheet("constantes")

    # --- aba 'constantes' (visual) ---
    visuais = [
        ("Observação", "A planilha usa números literais pedidos nas fórmulas."),
        ("a (m)", 6378137.0),
        ("f", 1/298.257223563),
        ("coef. lat. geocêntrica", 0.993305619977094),
        ("coef. gravidade elipsóide s^2", "0.0052790414, 0.0000232718, 0.0000001262, 0.0000000007"),
        ("coef. redução normal à altura", "0.00335281068118, 0.00344978600308"),
        ("correção atmosférica", "0.874, 9.9e-5, 3.5625e-9"),
        ("zero-tide (superfície)", "30.4, 91.2")
    ]
    ws_cst.write(0, 0, "Parâmetro")
    ws_cst.write(0, 1, "Valor")
    for i, (k, v) in enumerate(visuais, start=1):
        ws_cst.write(i, 0, k)
        if isinstance(v, (int, float)):
            ws_cst.write_number(i, 1, float(v))
        else:
            ws_cst.write(i, 1, str(v))

    # --- cabeçalhos da aba 'resultado' ---
    headers = [
        "ID",
        "Geodetic Lat (Tide Free)",
        "Geodetic Lon (Tide Free)",
        "h (Tide Free)",
        "g (mean tide)",
        "N (Zero Tide)",
        "Geocentric Lat (Tide Free)",            # G
        "h (Mean Tide)",                         # H
        "H (Zero Tide)",                         # I
        "H (Mean Tide)",                         # J
        "Atm Correction (Mean Tide)",            # K
        "g with atm correction (Mean Tide)",     # L
        "Normal gravity on the ellipsoid (Tide Free).",  # M
        "Normal Gravity on the surface (Mean Tide)",     # N
        "Normal gravity on the surface (Zero Tide).",    # O
        "g with atm correction (Zero Tide)",     # P  <--- NOVA
        "Gravity disturbance (Zero Tide).",      # Q  = P - O
        "fonte"                                  # R  <--- cópia
    ]
    for j, h in enumerate(headers):
        ws_res.write(0, j, h)

    # --- dados de entrada necessários ---
    cols_in = [
        "ID",
        "Geodetic Lat (Tide Free)",
        "Geodetic Lon (Tide Free)",
        "h (Tide Free)",
        "g (mean tide)",
        "N (Zero Tide)",
        "fonte"
    ]
    missing = [c for c in cols_in if c not in out_df.columns]
    if missing:
        raise KeyError(f"Faltam colunas na saída: {missing}")

    # escreve valores base (A..F) e a fonte (R)
    n = len(out_df)
    for i in range(n):
        r = i + 1
        ws_res.write(r, 0, out_df.iloc[i]["ID"])
        ws_res.write_number(r, 1, float(out_df.iloc[i]["Geodetic Lat (Tide Free)"]))
        ws_res.write_number(r, 2, float(out_df.iloc[i]["Geodetic Lon (Tide Free)"]))
        ws_res.write_number(r, 3, float(out_df.iloc[i]["h (Tide Free)"]))
        ws_res.write_number(r, 4, float(out_df.iloc[i]["g (mean tide)"]))
        ws_res.write_number(r, 5, float(out_df.iloc[i]["N (Zero Tide)"]))
        # coluna R (fonte)
        ws_res.write(r, 17, str(out_df.iloc[i]["fonte"]) if pd.notna(out_df.iloc[i]["fonte"]) else "")

    # --- fórmulas (usando funções em inglês: SIN/ATAN/TAN/PI) ---
    for i in range(2, n + 2):
        lat = f"$B${i}"  # Geodetic Lat (deg)
        h   = f"$D${i}"  # h tide free
        N   = f"$F${i}"  # N

        # G: geocentric latitude
        ws_res.write_formula(i-1, 6,  f"=ATAN(0.993305619977094*TAN({lat}*PI()/180))*180/PI()")

        # H: h mean tide
        ws_res.write_formula(i-1, 7,  f"={h}-(1+0.3-0.6)*(-0.198*((3/2)*(SIN(G{i}*PI()/180)^2)-0.5))")

        # I: H zero tide
        ws_res.write_formula(i-1, 8,  f"=D{i}-{N}")

        # J: H mean tide (ajuste adicional)
        ws_res.write_formula(i-1, 9,  f"=I{i}-(-0.198*((3/2)*(SIN(G{i}*PI()/180)^2)-0.5))")

        # K: atmospheric correction
        ws_res.write_formula(i-1,10, f"=0.874-(9.9*10^-5)*J{i}+(3.5625*10^-9)*J{i}^2")

        # L: g with atm correction (Mean Tide)
        ws_res.write_formula(i-1,11, f"=E{i}+K{i}")

        # M: gamma0 (ellipsoid) in mGal
        ws_res.write_formula(i-1,12,
            f"=(9.7803267715*(1+0.0052790414*(SIN({lat}*PI()/180))^2+"
            f"0.0000232718*(SIN({lat}*PI()/180))^4+"
            f"0.0000001262*(SIN({lat}*PI()/180))^6+"
            f"0.0000000007*(SIN({lat}*PI()/180))^8))*10^5")

        # N: normal gravity on surface (Mean Tide)
        ws_res.write_formula(i-1,13,
            f"=M{i}*(1-2*H{i}*(1+0.00335281068118-2*0.00335281068118*(SIN({lat}*PI()/180))*(SIN({lat}*PI()/180))+0.00344978600308)/6378137+(3*H{i}^2)/6378137^2)")

        # O: normal gravity on surface (Zero Tide)
        ws_res.write_formula(i-1,14, f"=N{i}+(30.4-91.2*(SIN(G{i}*PI()/180))^2)*0.001")

        # P: g with atm correction (Zero Tide)  ← NOVA
        ws_res.write_formula(i-1,15, f"=L{i}+(30.4-91.2*(SIN(G{i}*PI()/180))^2)*0.001")

        # Q: gravity disturbance (Zero Tide)  = g_zero - gamma_zero
        ws_res.write_formula(i-1,16, f"=P{i}-O{i}")

    ws_res.freeze_panes(1, 1)
    ws_res.autofilter(0, 0, n, len(headers) - 1)


# =============================
# GUI
# =============================
def read_grav_file(path: str) -> pd.DataFrame:
    # Lê Excel/CSV e usa as 6 primeiras colunas na ordem:
    # ID, lat, lon, h, g, fonte
    try:
        df = pd.read_excel(path, header=0)
    except Exception:
        df = pd.read_csv(path, engine='python')

    df = df.iloc[:, :6].copy()
    df.columns = ['ID', 'lat_tide_free', 'lon_tide_free', 'h_tide_free', 'g_mean_tide', 'fonte']

    # normaliza tipos numéricos
    for c in ['lat_tide_free', 'lon_tide_free', 'h_tide_free', 'g_mean_tide']:
        df[c] = pd.to_numeric(df[c], errors='coerce')

    # normaliza longitude para −180..180 e garante Oeste negativo
    lon = df['lon_tide_free'].values.astype(float)
    lon = np.where(lon > 180.0, lon - 360.0, lon)
    if (lon > 0).mean() > 0.9:
        lon = -np.abs(lon)
    df['lon_tide_free'] = lon

    return df



class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Processador Gravimetria + Geóide")
        self.geometry("720x500")
        self.grav_path = tk.StringVar()
        self.grid_path = tk.StringVar()
        self.out_path  = tk.StringVar(value="GravityDisturbance_Output.xlsx")
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
        p=filedialog.askopenfilename()
        if p: self.grav_path.set(p)
    def browse_grid(self):
        p=filedialog.askopenfilename()
        if p: self.grid_path.set(p)
    def browse_out(self):
        p=filedialog.asksaveasfilename(defaultextension=".xlsx")
        if p: self.out_path.set(p)

    def run(self):
        try:
            grav_df = read_grav_file(self.grav_path.get())
            lons, lats, gridN = read_geoid_grid_lon_lat_N(self.grid_path.get())

            out_df = process_with_grid(grav_df, lons, lats, gridN)
            with pd.ExcelWriter(self.out_path.get(), engine='xlsxwriter') as writer:
                write_with_constants_sheet(out_df, writer)

            messagebox.showinfo("OK", f"Arquivo gerado:\n{self.out_path.get()}")
        except Exception as e:
            messagebox.showerror("Erro", str(e))


def main():
    app=App()
    app.mainloop()

if __name__=="__main__":
    main()
