import qrcode

# Email details
email_address = "gutsie@kuleuven.be"
subject = "Info Clinical Trial"
content = "Please, could you send me more information about your study?"

# Format the email URL
email_url = f"mailto:{email_address}?subject={subject}&body={content}"

# Generate the QR code
qr = qrcode.QRCode(version=1, error_correction=qrcode.constants.ERROR_CORRECT_M)
qr.add_data(email_url)
qr.make(fit=True)

# Create an image from the QR code
image = qr.make_image(fill_color="black", back_color="white")

# Save the QR code as an image file
qr_filename = "email_qr_code.png"
image.save(qr_filename)
