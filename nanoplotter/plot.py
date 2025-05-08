import os
from base64 import b64encode
from io import BytesIO
from urllib.parse import quote as urlquote
import sys
from kaleido.scopes.plotly import PlotlyScope
import logging
import plotly.io as pio


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
        """Return the base64 encoding of the figure file and insert in html image tag."""
        data_uri = b64encode(open(self.path, "rb").read()).decode("utf-8").replace("\n", "")
        return '<img src="data:image/png;base64,{0}">'.format(data_uri)

    def encode2(self):
        """Return the base64 encoding of the fig attribute and insert in html image tag."""
        buf = BytesIO()
        self.fig.savefig(buf, format="png", bbox_inches="tight", dpi=100)
        buf.seek(0)
        string = b64encode(buf.read())
        return '<img src="data:image/png;base64,{0}">'.format(urlquote(string))

    def save(self, settings):
        if not(self.only_report):
            if self.html:
                with open(self.path, "w") as html_out:
                    html_out.write(self.html)
                if not settings["no_static"]:
                    try:
                        for fmt in settings["format"]:
                            self.save_static2(fmt)
                    except (AttributeError, ValueError) as e:
                        p = os.path.splitext(self.path)[0] + ".png"
                        if os.path.exists(p):
                            os.remove(p)

                        logging.warning("No static plots are saved due to some kaleido problem:")
                        logging.warning(e)

            elif self.fig:
                # if settings["format"] is a list, save the figure in all formats
                if isinstance(settings["format"], list):
                    for fmt in settings["format"]:
                        self.fig.savefig(
                            fname=self.path + "." + fmt,
                            format=fmt,
                            bbox_inches="tight",
                        )
                else:
                    self.fig.savefig(fname=self.path, format=settings["format"], bbox_inches="tight")
            else:
                sys.exit("No method to save plot object: no html or fig defined.")

    def show(self):
        if self.fig:
            return self.fig.fig
        else:
            sys.stderr.write(".show not implemented for Plot instance without fig attribute!")

    def save_static(self, figformat):
        scope = PlotlyScope()
        with open(self.path.replace("html", figformat), "wb") as f:
            f.write(scope.transform(self.fig, format=figformat))
            logging.info(
                f"Saved {self.path.replace('.html', '')}  as {figformat} (or png for --legacy)"
            )
    
    
    def save_static2(self, figformat):
        json_fig = pio.read_json(self.path.replace("html", "json"))
        json_fig.write_image(self.path.replace("html", figformat), format=figformat)
        logging.info(
                f"Saved {self.path.replace('.html', '')}  as {figformat} (or png for --legacy)"
            )


    
    def save_json(self, json_path=None):
        """
        将 fig 对象保存为 JSON 文件。
        
        参数:
            json_path (str): JSON 文件的路径。如果未指定，则基于 self.path 构造路径。
        """
        if not self.fig:
            logging.error("无法保存 JSON 文件：fig 对象为空")
            return
        
        # 如果未指定 json_path，则基于 self.path 构造默认路径
        if json_path is None:
            json_path = self.path.replace(".html", ".json")  # 替换 .html 为 .json
        
        try:
            # 将 fig 对象导出为 JSON 数据
            json_data = self.fig.to_json()

            # 写入 JSON 文件
            with open(json_path, "w", encoding="utf-8") as f:
                f.write(json_data)

            logging.info(f"成功保存 JSON 文件: {json_path}")
        except Exception as e:
            logging.error(f"保存 JSON 文件失败: {e}")
