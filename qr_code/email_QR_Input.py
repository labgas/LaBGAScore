import qrcode

# Code for input
email_address = input("Enter email address: ")
subject = input("Enter subject: ")
content = input("Enter content: ")

# Command 
email_url = f"mailto:{email_address}?subject={subject}&body={content}"

# QR code generation
qr = qrcode.QRCode(version=1, error_correction=qrcode.constants.ERROR_CORRECT_M)
qr.add_data(email_url)
qr.make(fit=True)

# QR image
image = qr.make_image(fill_color="black", back_color="white")

# Save file
qr_filename = "input_QR_image.png"
image.save(qr_filename)
print(f"QR code generated and saved as '{qr_filename}'.")
