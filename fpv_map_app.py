import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import rasterio
import geopandas as gpd
from shapely.geometry import Point
import contextily as ctx
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.transforms import Affine2D


class FPVMapApp:
    def __init__(self, root):
        self.root = root
        self.root.title("FPV Heatmap Generator")
        self.root.geometry("1280x800")

        self.heightmap, self.transform, self.crs, self.bounds = None, None, None, None
        self.receiver_marker, self.receiver_px, self.receiver_py = None, None, None

        self.load_dem("elevation_map.tif")
        self.create_widgets()
        self.plot_initial_map()
        self.fig.canvas.mpl_connect('button_press_event', self.on_map_click)

    def load_dem(self, file_path):
        try:
            with rasterio.open(file_path) as src:
                self.heightmap = src.read(1)
                self.transform = src.transform
                self.crs = src.crs
                self.bounds = src.bounds
        except rasterio.errors.RasterioIOError:
            messagebox.showerror("Error", f"Could not find file: {file_path}")
            self.root.quit()

    def create_widgets(self):
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        self.fig = Figure(figsize=(10, 9), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=main_frame)
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        control_frame = ttk.Frame(main_frame, padding="10")
        control_frame.pack(side=tk.RIGHT, fill=tk.Y)

        ttk.Label(control_frame, text="Antenna Height (m):").pack(pady=(10, 0))
        self.antenna_height_var = tk.StringVar(value="2")
        ttk.Entry(control_frame, textvariable=self.antenna_height_var, width=10).pack()

        # שדה חדש: תדר
        ttk.Label(control_frame, text="Frequency (MHz):").pack(pady=(10, 0))
        self.freq_var = tk.StringVar(value="5800")
        ttk.Entry(control_frame, textvariable=self.freq_var, width=10).pack()

        # שדה חדש: Gain של אנטנת שידור (Tx)
        ttk.Label(control_frame, text="Transmitter Antenna Gain (dBi):").pack(pady=(10, 0))
        self.tx_gain_var = tk.StringVar(value="14")
        ttk.Entry(control_frame, textvariable=self.tx_gain_var, width=10).pack()

        # שדה חדש: Gain של אנטנת קליטה (Rx)
        ttk.Label(control_frame, text="Receiver Antenna Gain (dBi):").pack(pady=(10, 0))
        self.rx_gain_var = tk.StringVar(value="14")
        ttk.Entry(control_frame, textvariable=self.rx_gain_var, width=10).pack()

        ttk.Label(control_frame, text="Flight Altitude (m):").pack(pady=(10, 0))
        self.flight_altitude_var = tk.IntVar(value=50)
        self.altitude_label = ttk.Label(control_frame, text="50 m")
        self.altitude_label.pack()
        ttk.Scale(
            control_frame, from_=0, to=200, orient=tk.HORIZONTAL,
            variable=self.flight_altitude_var,
            command=lambda s: self.altitude_label.config(text=f"{int(float(s))} m")
        ).pack(fill=tk.X, pady=(0, 20))

        self.run_button = ttk.Button(control_frame, text="Calculate / Update", state=tk.DISABLED,
                                     command=self.run_simulation)
        self.run_button.pack(pady=20, ipady=5)

        ttk.Label(control_frame, text="Instructions:", font=("Arial", 10, "bold")).pack(pady=(20, 5))
        ttk.Label(control_frame, text="1. Click map to place receiver.\n"
                                      "2. Set antenna & flight altitude.\n"
                                      "3. Click 'Calculate / Update'.",
                  justify=tk.LEFT).pack()

    def plot_initial_map(self):
        self.ax.clear()
        self.ax.set_xlim(self.bounds.left, self.bounds.right)
        self.ax.set_ylim(self.bounds.bottom, self.bounds.top)
        ctx.add_basemap(self.ax, crs=self.crs.to_string(), source=ctx.providers.Esri.WorldImagery)
        self.ax.set_title("Click to Place Receiver")
        self.ax.set_aspect('equal', adjustable='box')
        self.fig.tight_layout()
        self.canvas.draw()
        if self.receiver_marker:
            self.receiver_marker.remove()
            self.receiver_marker = None

    def on_map_click(self, event):
        if event.inaxes != self.ax: return

        real_x, real_y = event.xdata, event.ydata
        self.receiver_px, self.receiver_py = [int(coord) for coord in ~self.transform * (real_x, real_y)]

        # If a calculation was already run, clear the old heatmap to place a new marker
        if self.run_button['text'] == "Calculate / Update" and not self.run_button['state'] == tk.DISABLED:
            self.plot_initial_map()

        if self.receiver_marker: self.receiver_marker.remove()

        self.receiver_marker = self.ax.plot(real_x, real_y, 'rP', markersize=10, markeredgewidth=2)[0]
        self.canvas.draw()

        self.run_button.config(state=tk.NORMAL)

    def run_simulation(self):
        self.run_button.config(state=tk.DISABLED, text="Calculating...")
        self.root.update_idletasks()

        flight_altitude = self.flight_altitude_var.get()
        antenna_height = float(self.antenna_height_var.get())
        tx_gain_db = float(self.tx_gain_var.get())
        rx_gain_db = float(self.rx_gain_var.get())
        freq_mhz = float(self.freq_var.get())

        tx_power_dbm = 25 + tx_gain_db + rx_gain_db  # gain של Tx וגם Rx מחושבים

        receiver_sensitivity_dbm = -85

        lambda_wave = 3e8 / (freq_mhz * 1e6)

        def knife_edge_loss(v):
            if v > -0.7: return 6.91 + 20 * np.log10(np.sqrt((v - 0.1) ** 2 + 1) + v - 0.1)
            return 0

        rssi_map = np.full_like(self.heightmap, -120.0, dtype=float)
        rows, cols = self.heightmap.shape
        receiver_ground_height = self.heightmap[self.receiver_py, self.receiver_px]
        receiver_abs_height = receiver_ground_height + antenna_height

        for y in range(rows):
            for x in range(cols):
                if x == self.receiver_px and y == self.receiver_py: continue
                drone_abs_height = self.heightmap[y, x] + flight_altitude
                num_steps = int(np.ceil(np.sqrt((x - self.receiver_px) ** 2 + (y - self.receiver_py) ** 2)))
                if num_steps < 2: num_steps = 2
                x_line, y_line = np.linspace(self.receiver_px, x, num_steps), np.linspace(self.receiver_py, y,
                                                                                          num_steps)
                profile_heights = self.heightmap[y_line.astype(int), x_line.astype(int)]
                pixel_width = self.transform.a
                profile_distances = np.sqrt(
                    (x_line - self.receiver_px) ** 2 + (y_line - self.receiver_py) ** 2) * pixel_width
                if not profile_distances.any(): continue
                total_distance_m = profile_distances[-1]
                if total_distance_m <= 0: continue
                fspl = 20 * np.log10(total_distance_m / 1000.0) + 20 * np.log10(freq_mhz) + 32.44
                line_of_sight = np.linspace(receiver_abs_height, drone_abs_height, len(profile_heights))
                clearance = line_of_sight - profile_heights
                diffraction_loss = 0
                if np.min(clearance) < 0:
                    h_index = np.argmin(clearance)
                    h = -clearance[h_index]
                    d1 = profile_distances[h_index]
                    d2 = total_distance_m - d1
                    if d1 > 0 and d2 > 0:
                        v = h * np.sqrt((2 / lambda_wave) * (total_distance_m / (d1 * d2)))
                        diffraction_loss = knife_edge_loss(v)
                rssi = tx_power_dbm - fspl - diffraction_loss
                rssi_map[y, x] = rssi

        self.plot_final_heatmap(rssi_map, receiver_sensitivity_dbm)
        self.run_button.config(state=tk.NORMAL, text="Calculate / Update")

    def plot_final_heatmap(self, rssi_map, threshold):
        # Clears the entire figure to prevent multiple colorbars and old plots
        self.fig.clear()
        self.ax = self.fig.add_subplot(111)

        # Set the geographic extent for the plot, this is crucial for imshow
        extent = [self.bounds.left, self.bounds.right, self.bounds.bottom, self.bounds.top]

        # 1. Draw the basemap on the bottom layer (zorder=0)
        self.ax.set_xlim(self.bounds.left, self.bounds.right)
        self.ax.set_ylim(self.bounds.bottom, self.bounds.top)
        ctx.add_basemap(self.ax, crs=self.crs.to_string(), source=ctx.providers.Esri.WorldImagery, zorder=0)

        # Make weak signals transparent so they don't block the map view
        display_map = rssi_map.copy()
        display_map[display_map <= threshold] = np.nan

        # 2. Draw the heatmap as a continuous IMAGE using imshow (zorder=10)
        im = self.ax.imshow(display_map, cmap='inferno', alpha=0.3, extent=extent, zorder=10,
                            vmin=threshold, vmax=25)
        # Manually add a colorbar for the image
        self.fig.colorbar(im, ax=self.ax, label="Signal Strength (RSSI in dBm)", shrink=0.7)

        # 3. Draw the CONTOUR LINES on top of the heatmap (zorder=20)
        rows, cols = self.heightmap.shape
        x_coords = np.linspace(self.bounds.left, self.bounds.right, cols)
        y_coords = np.linspace(self.bounds.top, self.bounds.bottom, rows)  # Note: y-axis is inverted for images
        xx, yy = np.meshgrid(x_coords, y_coords)

        levels = np.arange(np.nanmin(self.heightmap), np.nanmax(self.heightmap), 10)
        self.ax.contour(xx, yy, self.heightmap, levels=levels, colors='white',
                        linewidths=0.7, alpha=0.8, zorder=20)

        # 4. Draw the RECEIVER MARKER on the very top layer (zorder=30)
        receiver_real_x, receiver_real_y = self.transform * (self.receiver_px, self.receiver_py)
        self.ax.plot(receiver_real_x, receiver_real_y, 'rx', markersize=12, markeredgewidth=2, zorder=30)

        self.ax.set_title('FPV Signal Quality Heatmap')
        self.canvas.draw()
if __name__ == "__main__":
    root = tk.Tk()
    app = FPVMapApp(root)
    root.mainloop()