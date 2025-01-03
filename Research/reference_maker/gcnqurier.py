#%%
import os
import csv
from datetime import datetime, timedelta
import gcn
from gcn.voevent import parse
#%%
class GCNAlertProcessor:
    """
    A class to process and store GCN alerts.
    """

    def __init__(self, output_file="alerts.csv", retention_period=30):
        """
        Initialize the GCN Alert Processor.

        Parameters:
            output_file (str): Path to the output CSV file.
            retention_period (int): Retention period for alerts in days.
        """
        self.output_file = output_file
        self.retention_period = retention_period

        # Initialize the CSV file with headers if it doesn't exist
        if not os.path.exists(self.output_file):
            with open(self.output_file, "w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow(["Object Name", "RA", "Dec", "Discovery Date", "Received At"])

    def handle_alert(self, payload):
        """
        Handle the received GCN alert.

        Parameters:
            payload (bytes): The raw VOEvent XML payload.
        """
        voevent = parse(payload)
        object_name, ra, dec, discovery_date = self.extract_alert_info(voevent)
        received_at = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")

        if object_name and ra and dec and discovery_date:
            print(f"Received Alert - Name: {object_name}, RA: {ra}, Dec: {dec}, Discovery Date: {discovery_date}")

            # Append the alert to the CSV file
            with open(self.output_file, "a", newline="") as file:
                writer = csv.writer(file)
                writer.writerow([object_name, ra, dec, discovery_date, received_at])

            # Clean up old alerts
            self.cleanup_old_alerts()

    def extract_alert_info(self, voevent):
        """
        Extract the relevant information from the VOEvent.

        Parameters:
            voevent (VOEvent): Parsed VOEvent object.

        Returns:
            tuple: Object name, RA, Dec, and discovery date (or None if not found).
        """
        try:
            object_name = next(
                (param.attrib["value"] for param in voevent.What.Param if param.attrib["name"] == "GRB Name"),
                None,
            )
            ra = next(
                (param.attrib["value"] for param in voevent.What.Param if param.attrib["name"] == "RA"),
                None,
            )
            dec = next(
                (param.attrib["value"] for param in voevent.What.Param if param.attrib["name"] == "Dec"),
                None,
            )
            discovery_date = next(
                (param.attrib["value"] for param in voevent.What.Param if param.attrib["name"] == "Discovery Date"),
                None,
            )
            return object_name, ra, dec, discovery_date
        except Exception as e:
            print(f"Error extracting alert info: {e}")
            return None, None, None, None

    def cleanup_old_alerts(self):
        """
        Remove alerts older than the retention period from the CSV file.
        """
        cutoff_date = datetime.utcnow() - timedelta(days=self.retention_period)
        updated_rows = []

        with open(self.output_file, "r") as file:
            reader = csv.reader(file)
            headers = next(reader)
            for row in reader:
                alert_date = datetime.strptime(row[4], "%Y-%m-%d %H:%M:%S")  # Received At column
                if alert_date >= cutoff_date:
                    updated_rows.append(row)

        # Overwrite the file with updated rows
        with open(self.output_file, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(headers)
            writer.writerows(updated_rows)

if __name__ == "__main__":
    # Initialize the processor
    processor = GCNAlertProcessor(output_file="alerts.csv", retention_period=30)

    # Listen for alerts
    print("Listening for GCN alerts...")
    listen(handler=processor.handle_alert)

# %%
#!/usr/bin/env python
import gcn
#%%
# Define your custom handler here.
@gcn.include_notice_types(
    gcn.notice_types.FERMI_GBM_FLT_POS,  # Fermi GBM localization (flight)
    gcn.notice_types.FERMI_GBM_GND_POS,  # Fermi GBM localization (ground)
    gcn.notice_types.FERMI_GBM_FIN_POS)  # Fermi GBM localization (final)
def handler(payload, root):
    # Look up right ascension, declination, and error radius fields.
    pos2d = root.find('.//{*}Position2D')
    ra = float(pos2d.find('.//{*}C1').text)
    dec = float(pos2d.find('.//{*}C2').text)
    radius = float(pos2d.find('.//{*}Error2Radius').text)

    # Print.
    print('ra = {:g}, dec={:g}, radius={:g}'.format(ra, dec, radius))

# Listen for VOEvents until killed with Control-C.
gcn.listen(handler=handler)
# %%
