import os
from base64 import b64encode
from io import BytesIO
from urllib.parse import quote as urlquote
import sys
from kaleido.scopes.plotly import PlotlyScope
import logging


class Plot(object):
    """A Plot object is defined by a path to the output file and the title of the plot."""

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
        """Return the base64 encoding of the figure file and insert in html image tag."""
        data_uri = b64encode(open(self.path, 'rb').read()).decode('utf-8').replace('\n', '')
        return '<img src="data:image/png;base64,{0}">'.format(data_uri)

    def encode2(self):
        """Return the base64 encoding of the fig attribute and insert in html image tag."""
        buf = BytesIO()
        self.fig.savefig(buf, format='png', bbox_inches='tight', dpi=100)
        buf.seek(0)
        string = b64encode(buf.read())
        return '<img src="data:image/png;base64,{0}">'.format(urlquote(string))

    def save(self, figformat):
        if self.html:
            with open(self.path, 'w') as html_out:
                html_out.write(self.html)
            try:
                self.save_static(figformat)
            except (AttributeError, ValueError) as e:
                p = os.path.splitext(self.path)[0]+".png"
                if os.path.exists(p):
                    os.remove(p)

                logging.warning("No static plots are saved due to some kaleido problem:")
                logging.warning(e)

        elif self.fig:
            self.fig.savefig(
                fname=self.path,
                format=figformat,
                bbox_inches='tight')
        else:
            sys.exit("No method to save plot object: no html or fig defined.")

    def show(self):
        if self.fig:
            return self.fig.fig
        else:
            sys.stderr.write(".show not implemented for Plot instance without fig attribute!")

    def save_static(self, figformat):
        scope = PlotlyScope()
        with open(self.path.replace('html', figformat), "wb") as f:
            f.write(scope.transform(self.fig, format=figformat))
            logging.info(
                f"Saved {self.path.replace('.html', '')}  as {figformat} (or png for --legacy)")
