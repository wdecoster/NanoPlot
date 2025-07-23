import os
from base64 import b64encode
from io import BytesIO
from urllib.parse import quote as urlquote
import sys
import logging

# bring in kaleido and ensure Chrome is installed
import kaleido
# this will download a small headless Chrome build the first time you run it
kaleido.get_chrome_sync()

from kaleido import write_fig_sync


class Plot(object):
    """A Plot object is defined by a path to the output file and the title of the plot."""

    only_report = False
    
    def __init__(self, path, title):
        self.path = path
        self.title = title
        self.fig = None
        self.html = None

    def encode(self):
        if self.html:
            return self.html
        elif self.fig:
            return self.encode2()
        else:
            return self.encode1()

    def encode1(self):
        data_uri = b64encode(open(self.path, "rb").read()).decode("utf-8").replace("\n", "")
        return '<img src="data:image/png;base64,{0}">'.format(data_uri)

    def encode2(self):
        buf = BytesIO()
        self.fig.savefig(buf, format="png", bbox_inches="tight", dpi=100)
        buf.seek(0)
        string = b64encode(buf.read())
        return '<img src="data:image/png;base64,{0}">'.format(urlquote(string))

    def save(self, settings):
        if not self.only_report:
            if self.html:
                with open(self.path, "w") as html_out:
                    html_out.write(self.html)
                if not settings["no_static"]:
                    try:
                        for fmt in settings["format"]:
                            self.save_static(fmt)
                    except (AttributeError, ValueError) as e:
                        p = os.path.splitext(self.path)[0] + ".png"
                        if os.path.exists(p):
                            os.remove(p)
                        logging.warning("No static plots are saved due to some kaleido problem:")
                        logging.warning(e)

            elif self.fig:
                if isinstance(settings["format"], list):
                    for fmt in settings["format"]:
                        self.fig.savefig(
                            fname=self.path + "." + fmt,
                            format=fmt,
                            bbox_inches="tight",
                        )
                else:
                    self.fig.savefig(
                        fname=self.path,
                        format=settings["format"],
                        bbox_inches="tight",
                    )
            else:
                sys.exit("No method to save plot object: no html or fig defined.")

    def show(self):
        if self.fig:
            return self.fig.fig
        else:
            sys.stderr.write(".show not implemented for Plot instance without fig attribute!")

    def save_static(self, figformat):
        """
        Export a Plotly figure via Kaleido v1’s write_fig_sync.
        """
        output_path = self.path.replace(".html", f".{figformat}")
        try:
            write_fig_sync(self.fig, path=output_path)
            logging.info(f"Saved {output_path} as {figformat}")
        except Exception as e:
            logging.warning("No static plots are saved due to some kaleido problem:")
            logging.warning(e)

