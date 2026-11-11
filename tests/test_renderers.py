import buttle.renderers as btr
import buttle.models.molecules as btmm
import buttle.models.alignments as btma
from PIL import Image

def test_save_alignment_image_png(tmp_path):
    output = tmp_path / "alignment.png"

    image = btr.render_alignment_image(
        btma.Profile(name=["1","2"], alignment=[list("ACGT"), list("ACCT")], molecule_type="dna"),
        scale = 1,
        wrap = 300,
        start = 1,
        consensus = False,
        diff = False,
    )
    btr.save_alignment_image(
        image,
        suffix = ".png",
        filename=str(output),
    )

    with Image.open(output) as image:
        assert image.format == "PNG"
        assert image.width > 0
        assert image.height > 0
        image.load()


def test_save_alignment_image_tif(tmp_path):
    output = tmp_path / "alignment.tiff"

    image = btr.render_alignment_image(
        btma.Profile(name=["1","2"], alignment=[list("ACGT"), list("ACCT")], molecule_type="dna"),
        scale = 4,
        wrap = 300,
        start = 1,
        consensus = False,
        diff = False,
    )
    btr.save_alignment_image(
        image,
        suffix = ".tiff",
        filename=str(output),
    )

    with Image.open(output) as image:
        assert image.format == "TIFF"
        assert image.width > 0
        assert image.height > 0
        image.load()
