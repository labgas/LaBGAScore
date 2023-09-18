import qrcode

# URL var: enter link/url for the pdf
pdf_url = "https://drive.google.com/u/0/uc?id=1Uvs7RSEbtPUsUvyRcrEEizRhAgqhlElH&export=download"

# QR code instance
qr = qrcode.QRCode(
    version=1,  # QR code version (adjust as needed)
    error_correction=qrcode.constants.ERROR_CORRECT_L,  # Error correction level
    box_size=10,  # Size of each box in the QR code
    border=4,     # Border size
)

# Add the PDF URL as data to the QR code
qr.add_data(pdf_url)
qr.make(fit=True)

# Create a PIL image from the QR code data
img = qr.make_image(fill_color="black", back_color="white")

# Save the QR code image as a PNG file
img.save("QR_to_pdf.png")