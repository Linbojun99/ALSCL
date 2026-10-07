"""Fetch a pinned, checksum-verified KaTeX distribution for self-contained output."""
from pathlib import Path
import hashlib
import io
import ssl
import tarfile
import urllib.request

VERSION = '0.16.22'
SHA256 = 'e9e0d167db3175481cbadaff38e8d90b130f6a3ddb451a47e43c577fd511f365'

def prepare_katex(destination):
    destination=Path(destination)
    marker=destination/'.verified'
    if marker.exists() and marker.read_text().strip()==SHA256: return
    context=ssl.create_default_context()
    # Python.org macOS installs may not include the system root certificates.
    if Path('/etc/ssl/cert.pem').exists(): context.load_verify_locations('/etc/ssl/cert.pem')
    with urllib.request.urlopen(f'https://registry.npmjs.org/katex/-/katex-{VERSION}.tgz',context=context,timeout=60) as response:
        data=response.read()
    if hashlib.sha256(data).hexdigest()!=SHA256: raise ValueError('KaTeX archive checksum mismatch')
    with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as archive:
        for member in archive.getmembers():
            if not member.isfile():continue
            if member.name in ('package/dist/katex.min.css','package/dist/katex.min.js','package/LICENSE') or member.name.startswith('package/dist/fonts/'):
                relative=member.name.removeprefix('package/dist/').removeprefix('package/')
                path=destination/relative
                if not path.resolve().is_relative_to(destination.resolve()):raise ValueError('Invalid archive path')
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_bytes(archive.extractfile(member).read())
    marker.write_text(SHA256+'\n')
