設定と作業支援
================================

設定、負荷見積もり、密充填モデルの生成、runops adapter は numpy や matplotlib を使わずに import できます。
設定の操作手順は `設定を作成・検証する <https://nkzono99.github.io/BEACH/configuration.html>`_ を参照してください。

設定
--------

.. automodule:: beach.config
   :members:
   :imported-members:
   :undoc-members:

負荷見積もり
------------------------

.. automodule:: beach.tooling.workload
   :members:

密充填モデル
------------------------

.. automodule:: beach.tooling.closepack_config
   :members:

runops adapter
--------------

.. autoclass:: beach.tooling.runops_adapter.Adapter
   :members:
